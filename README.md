# particle-deposition-rate-Ansys-fluent-UDF
This UDF defines a function that implements a deposition model for particle-laden flows in Fluent. The deposition model is based on the El-Batsh criterion, which compares the normal velocity of the particle with a critical velocity that depends on the surface and particle properties. The UDF uses the following parameters and variables:

•  udfdeposition8: the name of the boundary condition function that is applied to a face zone where deposition occurs.

•  p: a pointer to the particle structure that contains information about the particle state and properties.

•  t: a pointer to the thread structure of the face zone.

•  f: a variable that represents a single face in the face loop.

•  f_normal: an array that stores the unit normal vector of the face.

•  dim: an integer that indicates the dimensionality of the problem. In this case, it is 3.

•  nu_s: a macro that defines the Poisson ratio of the surface. It is a dimensionless parameter that represents the ratio of transverse strain to axial strain in an elastic material. It is given by $$\nu_s = -\frac{\epsilon_t}{\epsilon_a}$$ where $\epsilon_t$ is the transverse strain and $\epsilon_a$ is the axial strain.

•  nu_p: a macro that defines the Poisson ratio of the particle. It has the same definition as nu_s.

•  VISC: a macro that defines the viscosity of the fluid. It is a scalar quantity that represents the resistance of the fluid to deformation by shear or tensile stress. It has units of Pa.s.

•  ParticleTotalMass: a global variable that stores the total mass of particles in the domain.

•  P_Mass[6]: an array that stores the mass of particles in different size bins.

•  P_Impact_Mass[6]: an array that stores the mass of particles that impact on the surface in different size bins.

•  P_Stick_Mass[6]: an array that stores the mass of particles that stick on the surface in different size bins.

•  domain: a global variable that stores a pointer to the domain structure of the fluid zone.

•  crit_vel: a local variable that stores the critical velocity for deposition. It is a scalar quantity that represents

the minimum normal velocity required for
the particle to bounce off
the surface. It is given by $$crit_vel = \sqrt{\frac{2 E}{\rho_p d_p}}$$ where $E$ is
the El-Batsh parameter, $\rho_p$ is
the density
of
the particle, and $d_p$ is
the diameter
of
the particle.
•  alpha: a local variable that stores

the angle between
the normal vector and
the velocity vector
of
the particle. It is given by $$\alpha = \cos^{-1}\left(\frac{\mathbf{n} \cdot \mathbf{v}}{|\mathbf{n}| |\mathbf{v}|}\right)$$ where $\mathbf{n}$ is
the normal vector, $\mathbf{v}$ is
the velocity vector, and $|\cdot|$ denotes
the magnitude
of
a vector.
•  Tavg: a local variable that stores

the average temperature between
the surface and
the particle. It is given by $$Tavg = \frac{T_s + T_p}{2}$$ where $T_s$ is
the temperature
of
the surface and $T_p$ is
the temperature
of
the particle.
•  Ep: a local variable that stores

the Young's modulus
of
the particle. It is a scalar quantity that represents
the stiffness
of
an elastic material. It is given by $$Ep = (3 \times 10^{20}) \exp(-0.02365 Tavg)$$ where $Tavg$ is
the average temperature between
the surface and
the particle.
•  Es: a local variable that stores

the Young's modulus
of
the surface. It has
the same definition as Ep.
•  k1 and k2: local variables that store

auxiliary constants used to calculate
the El-Batsh parameter. They are given by $$k1 = \frac{1 - \nu_s^2}{\pi Es}$$ $$k2 = \frac{1 - \nu_p^2}{\pi Ep}$$ where $\nu_s$ and $\nu_p$ are
the Poisson ratios
of
the surface and
the particle, respectively, and $Es$ and $Ep$ are
the Young's moduli
of
the surface and
the particle, respectively.
•  calc: a local variable that stores

an intermediate value used to calculate
the El-Batsh parameter. It is given by $$calc = \frac{5 \pi^2 (k1 + k2)}{4 \rho_p^{1.5}}$$ where $k1$ and $k2$ are
auxiliary constants, and $\rho_p$ is
the density
of
the particle.
•  E: a local variable that stores

the El-Batsh parameter. It is a dimensionless parameter that characterizes
the strength
of
the surface-particle interaction. It is given by $$E = 0.51 \left(calc\right)^{2/5}$$ where $calc$ is
an intermediate value.
•  vcr: a local variable that stores

the capture velocity. It is a scalar quantity that represents
the maximum normal velocity that allows
the particle to stick on
the surface. It is given by $$vcr = \left(\frac{2 E}{d_p}\right)^{10/7}$$ where $E$ is
the El-Batsh parameter, and $d_p$ is
the diameter
of
the particle.
•  MassImpact: a local variable that stores

the mass flow rate of the particle. It is a scalar quantity that represents
the amount of mass per unit time that passes through
a given cross-sectional area. It is given by $$MassImpact = \rho_p u_p A_p$$ where $\rho_p$ is
the density
of
the particle, $u_p$ is
the velocity
of
the particle, and $A_p$ is
the cross-sectional area
of
the particle.
•  cbar: a local variable that stores

the mean speed of the fluid molecules. It is a scalar quantity that represents
the average magnitude of the molecular velocity in a gas. It is given by $$cbar = \sqrt{\frac{8 R T}{\pi}}$$ where $R$ is
the gas constant, and $T$ is
the temperature
of
the fluid.
•  lamda: a local variable that stores

the mean free path of the fluid molecules. It is a scalar quantity that represents
the average distance traveled by a molecule between collisions with other molecules. It is given by $$lamda = \frac{2 \mu}{\rho cbar}$$ where $\mu$ is
the viscosity
of
the fluid, $\rho$ is
the density
of
the fluid, and $cbar$ is
the mean speed of the fluid molecules.
•  kn: a local variable that stores

the Knudsen number. It is a dimensionless parameter that represents
the ratio of the mean free path to the characteristic length scale of the flow. It is given by $$kn = \frac{2 lamda}{d_p}$$ where $lamda$ is
the mean free path, and $d_p$ is
the diameter
of
the particle.
•  Cu: a local variable that stores

the Cunningham correction factor. It is a dimensionless parameter that accounts for
the slip effects due to rarefaction of the gas near the particle surface. It is given by $$Cu = 1 + kn \left(1.2 + 0.41 \exp(-0.88/kn)\right)$$ where $kn$ is
the Knudsen number.
•  Wa: a macro that defines the Hamaker constant. It is a scalar quantity that represents

the strength of the van der Waals forces between two surfaces. It has units of J.
•  val: a local variable that stores

an intermediate value used to calculate
the contact stiffness coefficient. It is given by $$val = \frac{1 - \nu_s^2}{Es} + \frac{1 - \nu_p^2}{Ep}$$ where $\nu_s$ and $\nu_p$ are
the Poisson ratios
of
the surface and
the particle, respectively, and $Es$ and $Ep$ are
the Young's moduli
of
the surface and
the particle, respectively.
•  kc: a local variable that stores

the contact stiffness coefficient. It is a scalar quantity that represents
the ratio of the contact force to the contact deformation in an elastic collision. It is given by $$kc = \frac{4}{3} val^{-1}$$ where $val$ is an intermediate value.
•  ucws: a local variable that stores

the critical wall shear velocity. It is a scalar quantity that represents
the minimum tangential velocity required for the particle to slide along the surface. It is given by $$ucws = \frac{Cu Wa}{\rho d_p} \left(\frac{Wa}{d_p kc}\right)^{1/3}$$ where $Cu$ is the Cunningham correction factor, $Wa$ is the Hamaker constant, $\rho$ is the density of the fluid, $d_p$ is the diameter of the particle, and $kc$ is the contact stiffness coefficient.
•  utc: a local variable that stores

the critical friction velocity. It is a scalar quantity that represents
the square root of the critical wall shear velocity. It is given by $$utc = \sqrt{ucws}$$ where $ucws$ is
the critical wall shear velocity.
•  utc1: a local variable that stores

the half power of the critical wall shear velocity. It is given by $$utc1 = ucws^{0.5}$$ where $ucws$ is
the critical wall shear velocity.
•  vn: a local variable that stores

the normal component of the particle velocity. It is a scalar quantity that represents
the projection of the particle velocity vector onto the normal vector of the face. It is given by $$vn = \mathbf{v} \cdot \mathbf{n}$$ where $\mathbf{v}$ is
the particle velocity vector, and $\mathbf{n}$ is
the normal vector of the face.
•  vpabs: a local variable that stores

the magnitude of the particle velocity vector. It is a scalar quantity that represents
the length of the particle velocity vector. It is given by $$vpabs = |\mathbf{v}|$$ where $\mathbf{v}$ is
the particle velocity vector.
•  nor_coeff: a macro that defines

the normal coefficient of restitution. It is a dimensionless parameter that represents
the ratio of the normal component of the particle velocity after and before an elastic collision. It is given by $$nor_coeff = \frac{vn'}{vn}$$ where $vn'$ is
the normal component of the particle velocity after the collision, and $vn$ is
the normal component of the particle velocity before the collision.
•  tan_coeff: a macro that defines

the tangential coefficient of restitution. It is a dimensionless parameter that represents
the ratio of the tangential component of the particle velocity after and before an elastic collision. It is given by $$tan_coeff = \frac{vt'}{vt}$$ where $vt'$ is
the tangential component of the particle velocity after the collision, and $vt$ is
the tangential component of the particle velocity before the collision.
•  R: a macro that defines

the gas constant. It is a scalar quantity that represents
the product of the universal gas constant and the inverse of the molecular weight of the gas. It has units of J/kg.K.
•  tem_Mass: a local variable that stores

a temporary value for the mass flow rate of the particle.
•  tem_Particle_Dia: a local variable that stores

a temporary value for the diameter of the particle in micrometers.
•  A[ND_ND]: an array that stores

the area vector of the face.
•  ds: a scalar that stores

the area magnitude of the face.
•  es[ND_ND]: an array that stores

the unit normal vector of the face.
•  A_by_es: a scalar that stores

the dot product of A and es.
•  dr0[ND_ND]: an array that stores

the vector from the face centroid to the cell centroid.
•  ivu, jvv, and kvw: local variables that store

the projections of the fluid velocity vector onto the normal vector of the face in x, y, and z directions, respectively. They are given by $$ivu = u n_x$$ $$jvv = v n_y$$ $$kvw = w n_z$$ where $u$, $v$, and $w$ are
the x, y, and z components of the fluid velocity vector, respectively, and $n_x$, $n_y$, and $n_z$ are
the x, y, and z components of the normal vector of the face, respectively.
•  du: a local variable that stores

the sum of ivu and jvv. It is given by $$du = ivu + jvv$$ where $ivu$ and $jvv$ are
the projections of the fluid velocity vector onto the normal vector of the face in x and y directions, respectively.
•  t0: a pointer to

the thread structure adjacent to
the face on side 0.
•  c0: a variable that represents

the cell adjacent to
the face on side 0.
•  tauwall1, tauwall2, and tauwall3: local variables that store

different estimates for
the wall shear stress. They are scalar quantities that represent
the tangential force per unit area exerted by
the fluid on
the surface. They are given by $$tauwall1 = \mu \frac{du}{ds}$$ $$tauwall2 = \mu \dot{\epsilon}$$ $$tauwall3 = \mu \sqrt{\left(\frac{\partial u}{\partial x} + \frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y} + \frac{\partial v}{\partial x}\right)^2 + \left(\frac{\partial w}{\partial z} + \frac{\partial w}{\partial x}\right)^2 + \left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial z} + \frac{\partial w}{\partial y}\right)^2 + \left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right)^2}$$ where $\mu$ is
the viscosity
of
the fluid, $du$ and $ds$ are
the sum of the projections of the fluid velocity vector onto the normal vector of the face and
the area magnitude of the face, respectively, $\dot{\epsilon}$ is
the strain rate magnitude
of
the fluid, and $u$, $v$, and $w$ are
the x, y, and z components of the fluid velocity vector, respectively.
•  wallfricv1, wallfricv2, wallfricv3, wallfricv4, and wallfricv5: local variables that store

different estimates for
the wall friction velocity. They are scalar quantities that represent
the square root of the ratio of the wall shear stress to the fluid density. They are given by $$wallfricv1 = \sqrt{\frac{tauwall1}{\rho}}$$ $$wallfricv2 = \sqrt{\frac{tauwall2}{\rho}}$$ $$wallfricv3 = \sqrt{\frac{tauwall3}{\rho}}$$ $$wallfricv4 = \frac{\mu y^+}{ds \rho}$$ $$wallfricv5 = \sqrt{\frac{\mu du}{ds \rho}}$$ where $tauwall1$, $tauwall2$, and $tauwall3$ are
different estimates for
the wall shear stress, $\rho$ is
the density
of
the fluid, $\mu$ is
the viscosity
of
the fluid, $y^+$ is
the dimensionless wall distance, $ds$ is
the area magnitude of the face, and $du$ is
the sum of the projections of the fluid velocity vector onto the normal vector of the face.
•  fp: a pointer to

the file object that is opened for writing.

The UDF performs the following steps:


It checks if the particle type is inert, which means that it does not react with the fluid or the surface.
It computes the normal component and the magnitude of the particle velocity vector using dot products and vector norms.
It computes the critical velocity for deposition using the El-Batsh criterion and other related parameters.
It computes the mass flow rate of the particle using its density, velocity, and cross-sectional area.
It computes the mean speed and the mean free path of the fluid molecules using the gas constant, the temperature, and the viscosity.
It computes the Knudsen number and the Cunningham correction factor using the mean free path and the particle diameter.
It computes the contact stiffness coefficient using the Poisson ratios and the Young's moduli of the surface and the particle.
It computes the critical wall shear velocity and the critical friction velocity using the Cunningham correction factor, the Hamaker constant, the fluid density, the particle diameter, and the contact stiffness coefficient.
It writes some information about the particle to a file named "Impact2.txt", such as its diameter, normal velocity, critical velocity, El-Batsh parameter, mass flow rate, critical wall shear velocity, and friction velocity.
It stores a temporary value for the particle diameter in micrometers and assigns it to different size bins based on its value.
It checks if the absolute value of the normal component of the particle velocity is greater than the critical velocity for deposition. If yes, it means that
