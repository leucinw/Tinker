# 
AMOEBA+ potential without charge flux implementation is fully supported in the main [Tinker release](https://github.com/TinkerTools/Tinker/tree/release).

This branch includes the preliminary implementation of charge flux into AMOEBA+ model.

We are currently merging charge flux implementation into the latest release.

Date: 11/22/2019


# AMOEBA+CF

  AMOEBA+ model was implemented based on the Tinker 8.2. Here are new features comparing to AMOEBA:

## vdW term
  A new combining rule (W-H) for epsilon was added (kvdw.f) 

## charge transfer (CT) term 
  A non-bonded interaction with user defined exclusion rules (1-2, 1-3, etc). 
  This is implemented as an individual energy term in the same manner as AMOEBA vdW term. 
  Except that CT does not use induction factor for hydrogen (i.e. induction factor = 1.0).
  (ect0.f, ect1.f and ect3.f)

## polarization
  AMOEBA+ uses a new damping scheme for the direct induction, while keeps the Thole damping for the mutual part.(induce.f, epolar.f, epolar1.f, epolar3.f)

## charge penetration (CP)
  We use slightly different damping functions than HIPPO model. (empole.f, empole1.f, empole3.f) 

## charge flux (CF) 

	!! Supported in AMOEBA+CF, not in AMOEBA+

  First I looped over bond and angle to accumulate charge fluxes(chrgflux.f). 
  Then I added the accumulated charges to the force field defined charges to get the updated charges(rotpol.f). 
  Potential is accumulated in the empole1.f and epolar1.f, and then a chain rule term due to charge flux is calculated by calling echrgflux1.f. 

## ind&inp induced dipole
  Currently in AMOEBA we have two sets of induced dipoles. In AMOEBA+ CPU code, I have merged them and use only one set of dipole. We used two sets of scaling factors (polar-1x-intra and polar-1x-inter).

## Reference

1. Liu, C.; Piquemal, J.-P.; Ren, P., Implementation of Geometry Dependent Charge Flux into AMOEBA+ Potential.  *submitted*
1. Liu, C.; Piquemal, J.-P.; Ren, P., AMOEBA+ Classical Potential for Modeling Molecular Interactions. *J. Chem. Theory Comput.* **2019**, 15 (7), 4122-4139(__Vdw, CT, CP, Pol__)
1. Liu, C.; Qi, R.; Wang, Q.; Piquemal, J.-P.; Ren, P., Capturing Many-body Interactions with Classical Dipole Induction Models. *J. Chem. Theory Comput.*, **2017**, 13 (6), 2751-2761 (__Pol__)
1. Rackers, J. A.; Wang, Q.; Liu, C.; Piquemal, J.-P.; Ren, P.; Ponder, J. W., An optimized charge penetration model for use with the AMOEBA force field. *Phys. Chem. Chem. Phys.* **2016**, 19, 276-291 (__CP__)
