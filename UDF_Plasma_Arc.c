#include "udf.h"
#include "sg_udms.h"
#include "mem.h"

#define epsion 8.854e-12;
#define w 0.5
#define Wk 28 
#define R 8.31451 
#define EpsionRelative 1 
//#define dp 1           
#define Kb 1.38e-23 
#define e 1.6e-19   
#define pi 3.141592653 

real sigma=0;//electrical conductivity
#define Bt C_UDMI(cell,thread,4) /* Magnetic field intensity (Tesla) */ 
#define Rad C_UDMI(cell,thread,5) /* Radiation power*/ 
#define Ohm C_UDMI(cell,thread,6) /* Ohm heating power*/ 
#define Lw C_UDMI(cell,thread,7)  /* Lorentz force work power*/ 
#define Jx C_UDMI(cell,thread,10) /* current density (A/m^2) x*/ 
#define Jy C_UDMI(cell,thread,11)/* current density (A/m^2) y*/
#define Jz C_UDMI(cell,thread,12)/* current density (A/m^2) z*/
#define J C_UDMI(cell,thread,13)/* current density (A/m^2)*/
#define Ax C_UDSI(cell,thread,Ax0) /* vector potential x*/ 
#define Ay C_UDSI(cell,thread,Ay1)/* vector potential y*/ 
//#define Az C_UDSI(cell,thread,Az2)/* vector potential z*/ 
#define Axy C_UDSI_G(cell,thread,Ax0)[1] /* vector potential gradient */ 
//#define Axz C_UDSI_G(cell,thread,Ax0)[2]
#define Ayx C_UDSI_G(cell,thread,Ay1)[0]
//#define Ayz C_UDSI_G(cell,thread,Ay1)[2]
//#define Azx C_UDSI_G(cell,thread,Az2)[0]
//#define Azy C_UDSI_G(cell,thread,Az2)[1]
#define Vx C_UDSI_G(cell,thread,V)[0] /* scalar potential gradient */ 
#define Vy C_UDSI_G(cell,thread,V)[1]/* y */ 
//#define Vz C_UDSI_G(cell,thread,V)[2]/* z */ 
#define Sig C_UDSI_DIFF(cell,thread,0) /* electric conductivity */ 
#define Mu 4*3.1415927e-7 /**/
#define P (C_P(cell,thread)+101325.0)/101325
#define M2 40
#define M1 2
#define e_charge   1.602e-19
#define kB  1.38e-23 
#define A 3.0e+4
#define phi_e 2.63
#define phi_c 4.52
#define m_e 9.1e-31
#define epsilon 0.35
#define alpha 5.6705e-8
#define Vi 15.68
#define Vd 20
real v=0.1;
DEFINE_INIT(my_init_func,d)
{
Thread *t;
/* Integrate dissipation. */
cell_t c;
thread_loop_c(t,d)
{
begin_c_loop(c,t)
 C_UDMI(c,t,0)=0;//Jx
 C_UDMI(c,t,1)=0;//Jy
 C_UDMI(c,t,2)=0;//Jz
 C_UDMI(c,t,3)=0;
 end_c_loop(c,t)
}

}
enum
{
	V,Ax0,Ay1,Az2,N_REQUIRED_UDS
};

DEFINE_SOURCE(en, cell, thread, dS, eqn)
{
  real source;
  Jx=-Sig*Vx;
  Jy=-Sig*Vy;
//  Jz=-Sig*Vz;
  Ohm=Sig*(pow(Vx,2)+pow(Vy,2));
  source=Ohm;
  dS[eqn]=0;
  C_UDMI(cell,thread,7)=source;
  return source*10;
}

real net(real x) 
{
	real y;
	if (x<2200)
		y=0;
	else if (x<3000)
		y= ((x - 2200) / (3000 - 2200)) * (116);
	else if (x<10000)
		y= 116 + ((x - 3000) / (10000 - 3000)) * (6659322493.1116 - 116);
	else if (x< 20000)
		y= 6659322493.1116 + ((x - 10000) / (20000 - 10000)) * (133240605237.293 - 6659322493.1116);
	else if (x< 30000)
		y= 133240605237.293 + ((x - 20000) / (30000 - 20000)) * (788312642696.99 - 133240605237.293);
	else
		y= 788312642696.99 + 1.4e8*(x - 30000);
	return y;

}

DEFINE_SOURCE(enegy_source, cell, thread, dS, eqn)
{
	real source;
	real Joule_effect;
	real radiation_loss;
	real T=C_T(cell, thread);
	  Jx=-Sig*Vx;
  Jy=-Sig*Vy;
//  Jz=-Sig*Vz;
	J=sqrt(pow(Jx,2)+pow(Jy,2));
	Lw = 2.5*kB*(Jx*C_T_G(cell, thread)[0] + Jy*C_T_G(cell, thread)[1]) / e_charge;
	Joule_effect = (pow(Jx, 2) + pow(Jy, 2)) / Sig;
	radiation_loss = 4. * M_PI * (/*PR*/net(T));
	source = Joule_effect + Lw - radiation_loss;
	if(source<0)source=0;
	dS[eqn] = 0;
	C_UDMI(cell,thread,10)=source;
	return source;

}

DEFINE_DIFFUSIVITY(uds0_diff,c,t,i)
{
	real T=C_T(c,t);

	if(T<=3000)
		sigma= 1;
	else if (T<=4300)
		sigma = 2;
	else if (T <= 5000)
		sigma = 12;
	else if (T <= 6000)
		sigma = 100;
	else if (T <= 7000)
		sigma = 125.24 + ((T - 6000) / (7000 - 6000)) * (552.6 - 125.24);
	else if (T <= 8000)
		sigma = 552.6 + ((T - 7000) / (8000 - 7000)) * (1309.22 - 552.6);
	else if (T <= 9000)
		sigma = 1309.22 + ((T - 8000) / (9000 - 8000)) * (2250 - 1309.22);
	else if (T <= 10000)
		sigma = 2250 + ((T - 9000) / (10000 - 9000)) * (3264 - 2250);
	else if (T <= 12000)
		sigma = 3264 + ((T - 10000) / (12000 - 10000)) * (5255 - 3264);
	else if (T <= 15000)
		sigma = 5255 + ((T - 12000) / (15000 - 12000)) * (7863 - 5255);
	else if (T <= 19000)
		sigma = 7863 + ((T - 15000) / (19000 - 15000)) * (11762 - 7863);
	else if (T <= 24000)
		sigma = 11762 + ((T - 19000) / (24000 - 19000)) * (15152 - 11762);
	else if (T <= 30000)
		sigma = 15152 + ((T - 19000) / (24000 - 19000)) * (17384 - 15152);
	else
		sigma= 17384 + 0.3 * (T-30000);
//	sigma=1.0;
return sigma;
}

DEFINE_SOURCE(x_velocity_source, cell, thread, dS, eqn)
{ 
     real source;
	 real x[ND_ND];
  	 real Bt1, Bt2;
   	 Bt1=C_UDSI_G(cell,thread,Ay1)[0]-C_UDSI_G(cell,thread,Ax0)[1];
//	 Bt2=C_UDSI_G(cell,thread,Ax0)[2]-C_UDSI_G(cell,thread,Ay1)[0];
	 source=Jy*Bt1;
	 C_CENTROID(x,cell,thread);

     dS[eqn]=0;

	 C_UDMI(cell, thread, 22)=source;

    return source;
}

DEFINE_SOURCE(y_velocity_source, cell, thread, dS, eqn)
{ 
     real source;
	 real x[ND_ND];
	 real Bt1, Bt2;
//   	 Bt1=C_UDSI_G(cell,thread,Az2)[1]-C_UDSI_G(cell,thread,Ay1)[2];
	 Bt2=C_UDSI_G(cell,thread,Ay1)[0]-C_UDSI_G(cell,thread,Ax0)[1];
	 source=-Jx*Bt2;
	 C_CENTROID(x,cell,thread);

     dS[eqn]=0;
	 	 C_UDMI(cell, thread, 23)=source;
    return source;
}

DEFINE_SOURCE(z_velocity_source, cell, thread, dS, eqn)
{ 
     real source;
	 real Bt1, Bt2;
   	 Bt1=C_UDSI_G(cell,thread,Ax0)[2]-C_UDSI_G(cell,thread,Az2)[0];
	 Bt2=C_UDSI_G(cell,thread,Az2)[1]-C_UDSI_G(cell,thread,Ay1)[2];
	 source=Jx*Bt1-Jy*Bt2;
     dS[eqn]=0;
	 if(source>=0)
		 source=0;
	 C_UDMI(cell, thread, 24)=source;
    return source;
}

DEFINE_SOURCE(heat_source,cell,thread,dS,eqn)
{
	real source;
	real T=C_T(cell,thread);
  real x[ND_ND];
  C_CENTROID(x,cell,thread);
	 if(x[0]<0)
{	 C_UDMI(cell,thread,0)=-C_UDSI_G(cell,thread,0)[0]*Sig;//Jx
 C_UDMI(cell,thread,1)=-C_UDSI_G(cell,thread,0)[1]*Sig;//Jy
  C_UDMI(cell,thread,2)=-C_UDSI_G(cell,thread,0)[2]*Sig;//Jx
	 }
	 else
	 {
		 C_UDMI(cell,thread,0)=C_UDSI_G(cell,thread,0)[0]*Sig;//Jx
 C_UDMI(cell,thread,1)=C_UDSI_G(cell,thread,0)[1]*Sig;//Jx
  C_UDMI(cell,thread,2)=C_UDSI_G(cell,thread,0)[2]*Sig;//Jx
	 }
 C_UDMI(cell,thread,3)=sqrt(C_UDMI(cell,thread,0)*C_UDMI(cell,thread,0)+C_UDMI(cell,thread,1)*C_UDMI(cell,thread,1));//Jx
	C_UDMI(cell,thread,4)=-C_UDSI_G(cell,thread,0)[0];
	C_UDMI(cell,thread,5)=-C_UDSI_G(cell,thread,0)[1];
	C_UDMI(cell,thread,6)=sqrt(C_UDMI(cell,thread,4)*C_UDMI(cell,thread,4)+C_UDMI(cell,thread,5)*C_UDMI(cell,thread,5));
source = C_UDMI(cell,thread,3)*C_UDMI(cell,thread,3)/Sig - 4. * 3.1415926 * (net(T));
dS[eqn] = 0;
if(sqrt(C_UDSI_G(cell,thread,0)[0]*C_UDSI_G(cell,thread,0)[0]+C_UDSI_G(cell,thread,0)[1]*C_UDSI_G(cell,thread,0)[1])>0)
C_UDMI(cell,thread,15)=1.167*C_R(cell,thread)/sqrt(C_UDSI_G(cell,thread,0)[0]*C_UDSI_G(cell,thread,0)[0]+C_UDSI_G(cell,thread,0)[1]*C_UDSI_G(cell,thread,0)[1]);
else
C_UDMI(cell,thread,15)=0;

C_UDMI(cell,thread,16)=sqrt(C_UDSI_G(cell,thread,0)[0]*C_UDSI_G(cell,thread,0)[0]+C_UDSI_G(cell,thread,0)[1]*C_UDSI_G(cell,thread,0)[1]);

C_UDMI(cell, thread, 24) = 4. * 3.1415926 * (net(T));
C_UDMI(cell,thread,25)=source;
return source;
}

DEFINE_SOURCE(Ax_source, cell, thread, dS, eqn)
{
  real source;
  real x[ND_ND];
  C_CENTROID(x,cell,thread);

  Jx=-Sig*Vx;

  source=Mu*Jx;
  dS[eqn]=0;
  return source;
}

DEFINE_SOURCE(Ay_source, cell, thread, dS, eqn)
{
  real source,r, x[ND_ND];
	C_CENTROID(x, cell, thread);

  Jy=-Sig*Vy;

  source=Mu*Jy;
  dS[eqn]=0;
  return source;
}

DEFINE_PROFILE(pressure_profile1,t,i)
{
real time; /* this will hold the position vector */
face_t f;
begin_f_loop(f,t)
{
time=CURRENT_TIME;
if(time<=0.004)
	F_PROFILE(f,t,i) = 0;
else
F_PROFILE(f,t,i) = -1.414*2685318*sin(100*time)*0.1;
}
end_f_loop(f,t)
}

DEFINE_DIFFUSIVITY(Amgeneral_diffusion, cell, thread, i)
{
	real x=C_T(cell,thread);
	real y;
	if (x <= 3500)
		y = 0.001;
	else if (x>3500&&x <= 6000)
		y=0.001+0.030728*x;
	else if (x>6000&&x <= 10000)
		y=76.82+0.728795*x;
	else if (x>10000&&x <= 15000)
		y=2990+0.99322*x;
	else if (x>15000&&x <= 20000)
		y=7956.1+0.66378*x;
	else if (x>20000&&x <= 24000)
	y=11275+0.115*x;
	else
		y=11735;
	return y;
}

/* Dependent on the itinerary curve, modifications are required */
DEFINE_CG_MOTION(piston,dt,vel,omega,time,dtime)
{
Thread *t;
face_t f;
real force, dv;
/* reset velocities */
NV_S(vel, =, 0.0);
NV_S(omega, =, 0.0);
if (!Data_Valid_P())
return;
/* get the thread pointer for which this motion is defined */
t = DT_THREAD(dt);

vel[0] = v;
}

DEFINE_PROFILE(pressure_profile,t,i)
{
real x[ND_ND]; /* this will hold the position vector */
real xx,yy;
real r;
face_t f;
begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
xx = x[0]-0.01*CURRENT_TIME;
yy = x[1]+0.005;

r=sqrt(xx*xx+yy*yy);
if(r<=0.005)
F_PROFILE(f,t,i) =2e7;
else
F_PROFILE(f,t,i) =0;
}
end_f_loop(f,t)
}

DEFINE_CG_MOTION(move, dt, vel, omega, time, dtime)
{
	NV_S(vel, =, 0.0);
	NV_S(omega, =, 0.0);
	

	/* Dependent on the itinerary curve, modifications are required */
	if (time <= 0.001425)
	{
		vel[0] = -10;
	}
	else if (time <= 0.002057)
	{
		vel[0] = -10;
	}
	else if (time <= 0.004003)
	{
		vel[0] = -10;
	}
	else if (time <= 0.00637)
	{
		vel[0] = -10;
	}
	else if (time <= 0.008157)
	{
		vel[0] = -10;
	}


	else if (time <= 0.009208)
	{
		vel[0] = -10;
	}
	else if (time <= 0.010205)
	{
		vel[0] = -10;
	}
	else if (time <= 0.0122)
	{
		vel[0] = -10;
	}
	else if (time <= 0.014195)
	{
		vel[0] = -10;
	}
	else if (time <= 0.016086)
	{
		vel[0] = -10;
	}


	else if (time <= 0.017925)
	{
		vel[0] = -10;
	}
	else if (time <= 0.018451)
	{
		vel[0] = -10;
	}
	else if (time <= 0.019924)
	{
		vel[0] = -10;
	}
	else if (time <= 0.021293)
	{
		vel[0] = -10;
	}
	else if (time <= 0.02282)
	{
		vel[0] = -10;
	}


	else if (time <= 0.025032)
	{
		vel[0] = -10;
	}
	else if (time <= 0.027244)
	{
		vel[0] = -10;
	}
	else if (time <= 0.03)
	{
		vel[0] = -10;
	}
	else
	{
		vel[0] = 0;
	}
	printf("\n");
	printf("\n x_velocity = %g \n", vel[0]);
}

DEFINE_PROFILE(jueyuan, thread, position)
{
	face_t f;
	cell_t c0;
	Thread* t0;
	real i0;
	real b = RP_Get_Real("flow-time");


/* Dependent on the current, modifications are required */
	
	i0 = (40000000) * 1.414 * sin(2 * 3.1415926 * 20 * b + 0.6 * 3.1415926);

	//Message(" time = %g \n",b);
   // Message(" CURRENT i0 = %g \n",i0);
	t0 = THREAD_T0(thread);
	begin_f_loop(f, thread)
	{
		c0 = F_C0(f, thread);
		C_UDMI(c0, t0, 27) = i0;
		F_PROFILE(f, thread, position) = i0;
	}
	end_f_loop(f, thread)
}



DEFINE_PROFILE(current_change,thread,position)
{
  face_t f;
  cell_t c0;
  Thread *t0;
  real i0;
  real b=RP_Get_Real("flow-time");

  /* Dependent on the current, modifications are required */

  i0 = (40000000) * 1.414 * sin(2 * 3.1415926 * 20 * b + 0.6 * 3.1415926);

 //Message(" time = %g \n",b);
// Message(" CURRENT i0 = %g \n",i0);
 t0=THREAD_T0(thread);
 begin_f_loop(f,thread)
    {
     c0=F_C0(f,thread);
     C_UDMI(c0,t0,27)=i0;
	 F_PROFILE(f, thread, position) = i0;
    }
 end_f_loop(f,thread)
}
