#include "udf.h"
#include "dpm.h"
#include "mem.h"
#include "sg.h"
#define nu_s 0.27 /*Poisson ratio of the surface*/
#define nu_p 0.27 /*Poisson ratio of the particle*/
#define VISC 1.7894e-05 /*Viscosity*/
/*#define RGAS (UNIVERSAL_GAS_CONSTANT/MW)*/
#define Tdatum 288.15
#define NUM_UDM 6 /* No. of user-defined memory locations8/
real ParticleTotalMass;
real P_Mass[6];
real P_Impact_Mass[6];
real P_Stick_Mass[6];
Domain *domain; *Get domain pointer and assign later to domain*/
/* Boundary condition macro for the deposition model*/
DEFINE_DPM_BC(udfdeposition8,p,t,f,f_normal,dim)
{
#if !RP_HOST /* Used only for 3D cases*/
real crit_vel,alpha,Tavg,Ep,vcr=0.,Es,k1,k2,calc,E,utau 1,dudz,wall_shear_stress;
real Cu,cbar,lamda,kn,val,kc,ucws,utc,h,del,ufr,ff,utau ,area,wall_shear_force,utaunew;
real vn=0.;
real utau2, utc1,vpabs=0.,MassImpact;
real nor_coeff=0.8;
real tan_coeff=0.4;
real R=287.;
real tem_Mass=0.;
real tem_Particle_Dia=0.;
real A[ND_ND],ds,es[ND_ND],A_by_es,dr0[ND_ND],ivu=0.,jvv=0.,kvw=0.,du=0.;
real tauwall1=0.,tauwall2=0.,tauwall3=0.,wallfricv1=0., wallfricv2=0.,wallfricv3=0.;
real yplus,wallfricv4=0.,wallfricv5=0.,uvel,vvel,wvel,d umag=0.;
FILE *fp;*
Thread *tcelll=P_CELL_THREAD(p);
cell_t c = P_CELL(p); /* Exact thread location of particle*/
/*cell_t c=RP_CELL(&p->cCell); */
/*Thread *tcell=RP_THREAD(&p->cCell);*/
cell_t c=P_CELL(p); /* Exact cell location of particle*/
Domain *d;
real normal[3];
int i,idim=dim;
real Wa = 0.039;
d=Get_Domain(1);
real x[ND_ND];
for (i=0.; i<idim; i++)
normal[i] = f_normal[i];
C_UDMI(c,tcell,0) += 1.0;
if(p->type==DPM_TYPE_INERT) /* Checks if inert particle is used*/
{
/*computing normal velocity*/
for(i=0.;i<idim;i++)
{
vn += p->state.V[i]*normal[i];
vpabs += pow(p->state.V[i],2.);
}
vpabs = pow(vpabs,.5);
/*computing critical velocity*/
Tavg = (F_T(f,t)+C_T(c,tcell))/2.;
/*Message("Avg temp is %g\n",Tavg);*/
/*Message("particle temp is %g\n",P_T(p));*/
/*Message("face temp from tcell is %g\n",F_T(f,tcell));*/
/*Message("normal face temp is %g\n",F_T(f,t));*/
/*NEW CORRELATION*/
Ep=(3.*(pow(10.,20.)))*exp(-0.02365*Tavg);
Es=Ep;
/*Message("Young's mod is %g\n",Ep);*/
k1 = (1.-(nu_s*nu_s))/(3.14*Es);
k2 = (1.-(nu_p*nu_p))/(3.14*Ep);
/*Message("particle density is %g\n",P_RHO(p));*/
/*Message("Ai particle density is %g\n",p->init_state.rho);*/
calc = (5.*3.14*3.14*(k1+k2))/(4.*(pow(P_RHO(p),1.5)));
E = 0.51*(pow(calc,(2./5.))); /* El-Batsh parameter*/
vcr = pow(((2.*E)/P_DIAM(p)),(10./7.));
/*Message("capture velocity is %g\n",vcr);*/
MassImpact=P_FLOW_RATE(p)*pow(10,9);
C_UDMI(c,tcell,1)+=MassImpact;
cbar = sqrt(((8.*C_RGAS(c,tcell)*C_T(c,tcell))/3.14));
lamda = (2.*C_MU_EFF(c,tcell))/(C_R(c,tcell)*cbar);
kn = (2.*lamda)/P_DIAM(p); /* Knudsen number*/
Cu = 1.+kn*(1.2+(0.41*exp(-0.88/kn))); /*Cunningham Correction Factor*/
Tavg = (F_T(f,t)+C_T(c,tcell))/2.; /*Avg of particle & surface temperature*/
/*Calculating El-Batsh parameter*/
val = ((1.-(pow(nu_s,2.)))/Es)+((1.-(pow(nu_p,2.)))/Ep);
kc = (4./3.)/val;
/*Calculating critical wall shear velocity*/
ucws = ((Cu*Wa)/(C_R(c,tcell)*P_DIAM(p)))*(pow((Wa/(P_DIAM(p)*kc)),(1./3.)));
utc = sqrt(ucws);
/* Message("wall shear vel is %g\n", utc);*/
utc1=pow(ucws,0.5);
/*Calculating wall friction velocity*/
/* Alternate Method 1*/
/*dudz = C_DUDZ(c,tcell);*/
/*utau1 = sqrt((VISC*dudz)/C_R(c,tcell));*/
/*Message("Friction vel by regular formula is %g\n",utau1);*/
yplus=F_STORAGE_R(f,t,SV_WALL_YPLUS_UTAU);
/*Alternate Method 2*/
ivu=C_U(c,tcell)*normal[0];
jvv=C_V(c,tcell)*normal[1];
kvw=C_W(c,tcell)*normal[2];
du =ivu+jvv;
Thread *t0=THREAD_T0(t);
cell_t c0=F_C0(f,t);
BOUNDARY_FACE_GEOMETRY(f,t,A,ds,es,A_by_es,dr0);
/* Message("wall distance is %g\n",ds);*/
tauwall1=C_MU_EFF(c,tcell)*(du/ds);
wallfricv1=sqrt(tauwall1/C_R(c,tcell));
/* Message("wallfrictionvelocity 1 is %g\n",wallfricv1);*/
if (yplus<11.25)
{
tauwall2=C_MU_EFF(c,tcell)*C_STRAIN_RATE_MAG(c,tce ll);
wallfricv2=sqrt(tauwall2/C_R(c,tcell));
}
else
{
wallfricv2 = (1./0.41)*log(yplus*9.);
}
/* Message("wallfrictionvelocity 2 is %g\n",wallfricv2);*/
/* Alternate Method 3*/
tauwall3=C_MU_EFF(c,tcell)*sqrt(C_DUDX(c,tcell)*(C _DUDX(c,tcell)+C_DUDX(c,tcell))+C_DUDY(c,tcell)*(C _DUDY(c,tcell)+C_DVDX(c,tcell))+C_DUDZ(c,tcell)*(C _D

UDZ(c,tcell)+C_DWDX(c,tcell))+C_DVDX(c,tcell)*(C_D VDX(c,tcell)+C_DUDY(c,tcell))+C_DVDY(c,tcell)*(C_D VDY(c,tcell)+C_DVDY(c,tcell))+C_DVDZ(c,tcell)*(C_D VDZ(c,tcell)+C_DWDY(c,tcell))+C_DWDX(c,tcell)*(C_D WDX(c,tcell)+C_DUDZ(c,tcell))+C_DWDY(c,tcell)*(C_D WDY(c,tcell)+C_DVDZ(c,tcell))+C_DWDZ(c,tcell)*(C_D WDZ(c,tcell)+C_DWDZ(c,tcell)));
wallfricv3=sqrt(tauwall3/C_R(c,tcell));
/*ACTUAL WALL FRICTION VELOCITY CALCULATION*/
/* Message("wall yplus is %g\n",yplus);*/
ds=C_WALL_DIST(c,tcell);
wallfricv4=(VISC*yplus)/(ds*C_R(c,tcell));
/* Alternate Method 5*/
uvel=C_U(c0,t0);
vvel=C_V(c0,t0);
wvel=C_W(c0,t0);
dumag=sqrt((uvel*uvel)+(vvel*vvel));
if(yplus>10)
wallfricv5=sqrt((C_MU_EFF(c,tcell)*dumag)/(C_R(c,tcell)*ds));
else
wallfricv5=sqrt((VISC*dumag)/(C_R(c,tcell)*ds));
/* Message("wallfrictionvelocity 5 is %g\n",wallfricv5);*/
fp=fopen("Impact2.txt","a");
fprintf(fp,"%6.3f %6.2f %6.2f %f %f %f %f\n", P_DIAM(p),vn,vcr,Ep,MassImpact,ucws,utaunew);
tem_Particle_Dia=P_DIAM(p)*pow(10,6);
if(tem_Particle_Dia<1)
P_Impact_Mass[0]+=P_FLOW_RATE(p)*pow(10,9);
else if(tem_Particle_Dia>1 && tem_Particle_Dia<3)
P_Impact_Mass[1]+=P_FLOW_RATE(p)*pow(10,9);
else if(tem_Particle_Dia>3 && tem_Particle_Dia<5)
P_Impact_Mass[2]+=P_FLOW_RATE(p)*pow(10,9);
else if(tem_Particle_Dia>5 && tem_Particle_Dia<7)
P_Impact_Mass[3]+=P_FLOW_RATE(p)*pow(10,9);
else if(tem_Particle_Dia>7 && tem_Particle_Dia<10)
P_Impact_Mass[4]+=P_FLOW_RATE(p)*pow(10,9);
else
P_Impact_Mass[5]+=P_FLOW_RATE(p)*pow(10,9);
/* PARTICLE DOES NOT STICK*/
if(abs(vn)>vcr)
{
C_UDMI(c,tcell,2) += 1.0;
C_UDMI(c,tcell,3) += P_FLOW_RATE(p)*pow(10,9);
for(i=0;i<idim;i++)
p->state.V[i]-=vn*normal[i];
/*Apply tangential coefficient of restitution. */
for (i=0; i<idim; i++)
p->state.V[i]*=tan_coeff;
/*Add reflected normal velocity. */
for (i=0; i<idim; i++)
p->state.V[i]-=nor_coeff*vn*normal[i];
/* Store new velocity in state0 of particle*/
for (i=0; i<idim; i++)
p->state0.V[i]=p->state.V[i];
return PATH_ACTIVE;
}
/* PARTICLE DEPOSITS*/
else
if(wallfricv4 < utc)
{
/*num of particles deposited or num of hits*/
C_UDMI(c,tcell,4) += 1.0;
/* mass of particles deposited*/
C_UDMI(c,tcell,5) += P_FLOW_RATE(p)*pow(10,9);
tem_Mass=P_FLOW_RATE(p)*pow(10,9);
if(tem_Particle_Dia<1.)
P_Stick_Mass[0]+=P_FLOW_RATE(p)*pow(10,9);
else if(tem_Particle_Dia>1 && tem_Particle_Dia<3)
P_Stick_Mass[1]+=P_FLOW_RATE(p)*pow(10,9);
else if(tem_Particle_Dia>3 && tem_Particle_Dia<5)
P_Stick_Mass[2]+=P_FLOW_RATE(p)*pow(10,9);
else if(tem_Particle_Dia>5 && tem_Particle_Dia<7)
P_Stick_Mass[3]+=P_FLOW_RATE(p)*pow(10,9);
else if(tem_Particle_Dia>7 && tem_Particle_Dia<10)
P_Stick_Mass[4]+=P_FLOW_RATE(p)*pow(10,9);
else
P_Stick_Mass[5]+=P_FLOW_RATE(p)*pow(10,9);
fp=fopen("Stick2.txt","a");
fprintf(fp,"%f %f %6.2f %f %6.2f %6.2f %f\n", tem_Particle_Dia,tem_Mass,vn,vcr,Ep,ucws,utaunew);
fclose(fp);
}
}
return PATH_ABORT;
#endif /* Only for 3D cases*/
}
/*Setting UDM locations and initializing them to 0*/
DEFINE_ON_DEMAND(reset_UDMsnew)
{
int i=0;
Thread *t;
Domain *d;
cell_t c;
face_t f;
d=Get_Domain(1);
Message("Setting UDMs \n");
thread_loop_c(t,d)
{
begin_c_loop(c,t)
{
for(i=0;i<6;i++)
C_UDMI(c,t,i)=0.0;
}
end_c_loop(c,t)
}
}