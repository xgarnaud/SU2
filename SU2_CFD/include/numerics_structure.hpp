/*!
 * \file numerics_structure.hpp
 * \brief Headers of the main subroutines for the dumerical definition of the problem.
 *        The subroutines and functions are in the <i>numerics_structure.cpp</i>,
 *        <i>numerics_convective.cpp</i>, <i>numerics_viscous.cpp</i>, and
 *        <i>numerics_source.cpp</i> files.
 * \author F. Palacios, T. Economon
 * \version 4.0.1 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include "../../Common/include/config_structure.hpp"
#include "numerics_machine_learning.hpp"
#include "numerics_machine_learning_turbulent.hpp"
#include "variable_structure.hpp"

using namespace std;

/*!
 * \class CNumerics
 * \brief Class for defining the numerical methods.
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CNumerics {
protected:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
  unsigned short nSpecies; 	/*!< \brief No of species present in plasma */
  su2double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
  su2double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
  su2double Gas_Constant;		 		/*!< \brief Gas constant. */
  su2double *Vector; /*!< \brief Auxiliary vector. */
  su2double *Enthalpy_formation;
  unsigned short nDiatomics, nMonatomics;
  su2double Prandtl_Lam;				/*!< \brief Laminar Prandtl's number. */
  su2double Prandtl_Turb;		/*!< \brief Turbulent Prandtl's number. */
  
public:
  
  su2double
  **Flux_Tensor,	/*!< \brief Flux tensor (used for viscous and inviscid purposes. */
  *Proj_Flux_Tensor;		/*!< \brief Flux tensor projected in a direction. */
  
  su2double
  **tau,		/*!< \brief Viscous stress tensor. */
  **delta;			/*!< \brief Identity matrix. */
  su2double **dVdU; /*!< \brief Transformation matrix from primitive variables, V, to conserved, U. */
  su2double
  *Diffusion_Coeff_i, /*!< \brief Species diffusion coefficients at point i. */
  *Diffusion_Coeff_j; /*!< \brief Species diffusion coefficients at point j. */
  su2double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
  Laminar_Viscosity_j,		/*!< \brief Laminar viscosity at point j. */
  Laminar_Viscosity_id,	/*!< \brief Variation of laminar viscosity at point i. */
  Laminar_Viscosity_jd;		/*!< \brief Variation of laminar viscosity at point j. */
  su2double Thermal_Conductivity_i, /*!< \brief Thermal conductivity at point i. */
  Thermal_Conductivity_j, /*!< \brief Thermal conductivity at point j. */
  Thermal_Conductivity_ve_i, /*!< \brief Thermal conductivity at point i. */
  Thermal_Conductivity_ve_j; /*!< \brief Thermal conductivity at point j. */
  su2double Cp_i, /*!< \brief Cp at point i. */
  Cp_j;         /*!< \brief Cp at point j. */
  su2double *Theta_v; /*!< \brief Characteristic vibrational temperature */
  su2double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
  Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
  su2double turb_ke_i,	/*!< \brief Turbulent kinetic energy at point i. */
  turb_ke_j;			/*!< \brief Turbulent kinetic energy at point j. */
  su2double Pressure_i,	/*!< \brief Pressure at point i. */
  Pressure_j;			/*!< \brief Pressure at point j. */
  su2double GravityForce_i,	/*!< \brief Gravity force at point i. */
  GravityForce_j;			/*!< \brief Gravity force at point j. */
  su2double Density_i,	/*!< \brief Density at point i. */
  Density_j;			/*!< \brief Density at point j. */
  su2double DensityInc_i,	/*!< \brief Incompressible density at point i. */
  DensityInc_j;			/*!< \brief Incompressible density at point j. */
  su2double BetaInc2_i,	/*!< \brief Beta incompressible at point i. */
  BetaInc2_j;			/*!< \brief Beta incompressible at point j. */
  su2double Lambda_i,	/*!< \brief Spectral radius at point i. */
  Lambda_j;			/*!< \brief Spectral radius at point j. */
  su2double LambdaComb_i,	/*!< \brief Spectral radius at point i. */
  LambdaComb_j;			/*!< \brief Spectral radius at point j. */
  su2double SoundSpeed_i,	/*!< \brief Sound speed at point i. */
  SoundSpeed_j;			/*!< \brief Sound speed at point j. */
  su2double Enthalpy_i,	/*!< \brief Enthalpy at point i. */
  Enthalpy_j;			/*!< \brief Enthalpy at point j. */
  su2double dist_i,	/*!< \brief Distance of point i to the nearest wall. */
  dist_j;			/*!< \brief Distance of point j to the nearest wall. */
  su2double Temp_i,	/*!< \brief Temperature at point i. */
  Temp_j;			/*!< \brief Temperature at point j. */
  su2double *Temp_tr_i, /*!< \brief Temperature transl-rot at point i. */
  *Temp_tr_j;/*!< \brief Temperature transl-rot at point j. */
  su2double *Temp_vib_i, /*!< \brief Temperature vibrational at point i. */
  *Temp_vib_j;/*!< \brief Temperature vibrational at point j. */
  su2double *Und_Lapl_i, /*!< \brief Undivided laplacians at point i. */
  *Und_Lapl_j;		/*!< \brief Undivided laplacians at point j. */
  su2double Sensor_i,	/*!< \brief Pressure sensor at point i. */
  Sensor_j;			/*!< \brief Pressure sensor at point j. */
  su2double *GridVel_i,	/*!< \brief Grid velocity at point i. */
  *GridVel_j;			/*!< \brief Grid velocity at point j. */
  su2double *U_i,		/*!< \brief Vector of conservative variables at point i. */
  *U_id,		/*!< \brief Vector of derivative of conservative variables at point i. */
  *UZeroOrder_i,  /*!< \brief Vector of conservative variables at point i without reconstruction. */
  *U_j,				/*!< \brief Vector of conservative variables at point j. */
  *UZeroOrder_j,  /*!< \brief Vector of conservative variables at point j without reconstruction. */
  *U_jd,				/*!< \brief Vector of derivative of conservative variables at point j. */
  *U_0,				/*!< \brief Vector of conservative variables at node 0. */
  *U_1,				/*!< \brief Vector of conservative variables at node 1. */
  *U_2,				/*!< \brief Vector of conservative variables at node 2. */
  *U_3;				/*!< \brief Vector of conservative variables at node 3. */
  su2double *V_i,		/*!< \brief Vector of primitive variables at point i. */
  *V_j;				/*!< \brief Vector of primitive variables at point j. */
  su2double *S_i,		/*!< \brief Vector of secondary variables at point i. */
  *S_j;				/*!< \brief Vector of secondary variables at point j. */
  su2double *Psi_i,		/*!< \brief Vector of adjoint variables at point i. */
  *Psi_j;				/*!< \brief Vector of adjoint variables at point j. */
  su2double *DeltaU_i,	/*!< \brief Vector of linearized variables at point i. */
  *DeltaU_j;			/*!< \brief Vector of linearized variables at point j. */
  su2double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
  *TurbVar_id,	/*!< \brief Vector of derivative of turbulent variables at point i. */
  *TurbVar_j,			/*!< \brief Vector of turbulent variables at point j. */
  *TurbVar_jd;	/*!< \brief Vector of derivative of turbulent variables at point j. */
  su2double *TransVar_i,	/*!< \brief Vector of turbulent variables at point i. */
  *TransVar_j;			/*!< \brief Vector of turbulent variables at point j. */
  su2double *LevelSetVar_i,	/*!< \brief Vector of turbulent variables at point i. */
  *LevelSetVar_j;			/*!< \brief Vector of turbulent variables at point j. */
  su2double *TurbPsi_i,	/*!< \brief Vector of adjoint turbulent variables at point i. */
  *TurbPsi_j;			/*!< \brief Vector of adjoint turbulent variables at point j. */
  su2double **ConsVar_Grad_i,	/*!< \brief Gradient of conservative variables at point i. */
  **ConsVar_Grad_j,			/*!< \brief Gradient of conservative variables at point j. */
  **ConsVar_Grad_0,			/*!< \brief Gradient of conservative variables at point 0. */
  **ConsVar_Grad_1,			/*!< \brief Gradient of conservative variables at point 1. */
  **ConsVar_Grad_2,			/*!< \brief Gradient of conservative variables at point 2. */
  **ConsVar_Grad_3,			/*!< \brief Gradient of conservative variables at point 3. */
  **ConsVar_Grad;				/*!< \brief Gradient of conservative variables which is a scalar. */
  su2double **PrimVar_Grad_i,	/*!< \brief Gradient of primitive variables at point i. */
  **PrimVar_Grad_j;			/*!< \brief Gradient of primitive variables at point j. */
  su2double *PrimVar_Lim_i,	/*!< \brief Limiter of primitive variables at point i. */
  *PrimVar_Lim_j;			/*!< \brief Limiter of primitive variables at point j. */
  su2double *PsiVar_Lim_i,		/*!< \brief Limiter of adjoint variables at point i. */
  *PsiVar_Lim_j;			/*!< \brief Limiter of adjoint variables at point j. */
  su2double **PsiVar_Grad_i,		/*!< \brief Gradient of adjoint variables at point i. */
  **PsiVar_Grad_j;			/*!< \brief Gradient of adjoint variables at point j. */
  su2double **TurbVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
  **TurbVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
  su2double **TransVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
  **TransVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
  su2double **LevelSetVar_Grad_i,	/*!< \brief Gradient of level set variables at point i. */
  **LevelSetVar_Grad_j;			/*!< \brief Gradient of level set variables at point j. */
  su2double **TurbPsi_Grad_i,	/*!< \brief Gradient of adjoint turbulent variables at point i. */
  **TurbPsi_Grad_j;			/*!< \brief Gradient of adjoint turbulent variables at point j. */
  su2double *AuxVar_Grad_i,		/*!< \brief Gradient of an auxiliary variable at point i. */
  *AuxVar_Grad_j;				/*!< \brief Gradient of an auxiliary variable at point i. */
  su2double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
  *Coord_j,			/*!< \brief Cartesians coordinates of point j. */
  *Coord_0,			/*!< \brief Cartesians coordinates of point 0 (Galerkin method, triangle). */
  *Coord_1,			/*!< \brief Cartesians coordinates of point 1 (Galerkin method, tetrahedra). */
  *Coord_2,			/*!< \brief Cartesians coordinates of point 2 (Galerkin method, triangle). */
  *Coord_3;			/*!< \brief Cartesians coordinates of point 3 (Galerkin method, tetrahedra). */
  unsigned short Neighbor_i,	/*!< \brief Number of neighbors of the point i. */
  Neighbor_j;					/*!< \brief Number of neighbors of the point j. */
  su2double *Normal,	/*!< \brief Normal vector, it norm is the area of the face. */
  *UnitNormal,		/*!< \brief Unitary normal vector. */
  *UnitNormald;		/*!< \brief derivatve of unitary normal vector. */
  su2double TimeStep,		/*!< \brief Time step useful in dual time method. */
  Area,				/*!< \brief Area of the face i-j. */
  Volume;				/*!< \brief Volume of the control volume around point i. */
  su2double Volume_n,	/*!< \brief Volume of the control volume at time n. */
  Volume_nM1,		/*!< \brief Volume of the control volume at time n-1. */
  Volume_nP1;		/*!< \brief Volume of the control volume at time n+1. */
  su2double *U_n,	/*!< \brief Vector of conservative variables at time n. */
  *U_nM1,		/*!< \brief Vector of conservative variables at time n-1. */
  *U_nP1;		/*!< \brief Vector of conservative variables at time n+1. */
  su2double vel2_inf; /*!< \brief value of the square of freestream speed. */
  su2double *WindGust_i,	/*!< \brief Wind gust at point i. */
  *WindGust_j;			/*!< \brief Wind gust at point j. */
  su2double *WindGustDer_i,	/*!< \brief Wind gust derivatives at point i. */
  *WindGustDer_j;			/*!< \brief Wind gust derivatives at point j. */
  su2double *Vorticity_i, *Vorticity_j;  /*!< \brief Vorticity. */
  su2double StrainMag_i, StrainMag_j;   /*!< \brief Strain rate magnitude. */
  
  su2double *l, *m;
  su2double *dPdU_i, *dPdU_j;
  su2double *dTdU_i, *dTdU_j;
  su2double *dTvedU_i, *dTvedU_j;
  su2double *Ys, **dFdYj, **dFdYi, *sumdFdYih, *sumdFdYjh, *sumdFdYieve, *sumdFdYjeve;
  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, VEL_INDEX, P_INDEX,
  RHO_INDEX, H_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;
  CVariable *var;
  
  /*!
   * \brief Constructor of the class.
   */
  CNumerics(void);
  
  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNumerics(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CNumerics(void);
  
  /*!
   * \brief Compute the determinant of a 3 by 3 matrix.
   * \param[in] val_matrix 3 by 3 matrix.
   * \result Determinant of the matrix
   */
  su2double Determinant_3x3(su2double A00, su2double A01, su2double A02,
                            su2double A10, su2double A11, su2double A12,
                            su2double A20, su2double A21, su2double A22);
  
  /*!
   * \brief Set the solution at different times.
   * \param[in] val_u_nM1 Conservative solution at time n-1.
   * \param[in] val_u_n Conservative solution at time n.
   * \param[in] val_u_nP1 Conservative solution at time n+1.
   */
  void SetPastSol(su2double *val_u_nM1, su2double *val_u_n, su2double *val_u_nP1);
  
  /*!
   * \brief Set the control volume at different times.
   * \param[in] val_volume_nM1 - Control volume at time n-1.
   * \param[in] val_volume_n - Control volume at time n.
   * \param[in] val_volume_nP1 - Control volume at time n+1.
   */
  void SetPastVolume(su2double val_volume_nM1, su2double val_volume_n, su2double val_volume_nP1);
  
  /*!
   * \brief Set the time step.
   * \param[in] val_timestep - Value of the time step.
   */
  void SetTimeStep(su2double val_timestep);
  
  /*!
   * \brief Get the Preconditioning Beta.
   * \return val_Beta - Value of the low Mach Preconditioner.
   */
  virtual su2double GetPrecond_Beta();
  
  /*!
   * \brief Set the freestream velocity square.
   * \param[in] SetVelocity2_Inf - Value of the square of the freestream velocity.
   */
  void SetVelocity2_Inf(su2double val_velocity2);
  
  /*!
   * \brief Set the value of the vorticity
   * \param[in] val_vorticity - Value of the vorticity.
   */
  void SetVorticity(su2double *val_vorticity_i, su2double *val_vorticity_j);
  
  /*!
   * \brief Set the value of the rate of strain magnitude.
   * \param[in] val_StrainMag_i - Value of the magnitude of rate of strain at point i.
   * \param[in] val_StrainMag_j - Value of the magnitude of rate of strain at point j.
   */
  void SetStrainMag(su2double val_strainmag_i, su2double val_strainmag_j);
  
  /*!
   * \brief Set the value of the conservative variables.
   * \param[in] val_u_i - Value of the conservative variable at point i.
   * \param[in] val_u_j - Value of the conservative variable at point j.
   */
  void SetConservative(su2double *val_u_i, su2double *val_u_j);
  
  /*!
   * \brief Set the value of the conservative variables withour reconstruction.
   * \param[in] val_u_i - Value of the conservative variable at point i.
   * \param[in] val_u_j - Value of the conservative variable at point j.
   */
  void SetConservative_ZeroOrder(su2double *val_u_i, su2double *val_u_j);
  
  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_v_i - Value of the primitive variable at point i.
   * \param[in] val_v_j - Value of the primitive variable at point j.
   */
  void SetPrimitive(su2double *val_v_i, su2double *val_v_j);
  
  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_v_i - Value of the primitive variable at point i.
   * \param[in] val_v_j - Value of the primitive variable at point j.
   */
  void SetSecondary(su2double *val_s_i, su2double *val_s_j);
  
  /*!
   * \brief Set the value of the conservative variables.
   * \param[in] val_u_0 - Value of the conservative variable at point 0.
   * \param[in] val_u_1 - Value of the conservative variable at point 1.
   * \param[in] val_u_2 - Value of the conservative variable at point 2.
   */
  void SetConservative(su2double *val_u_0, su2double *val_u_1, su2double *val_u_2);
  
  /*!
   * \brief Set the value of the conservative variables.
   * \param[in] val_u_0 - Value of the conservative variable at point 0.
   * \param[in] val_u_1 - Value of the conservative variable at point 1.
   * \param[in] val_u_2 - Value of the conservative variable at point 2.
   * \param[in] val_u_3 - Value of the conservative variable at point 3.
   */
  void SetConservative(su2double *val_u_0, su2double *val_u_1, su2double *val_u_2, su2double *val_u_3);
  
  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad_i - Gradient of the conservative variable at point i.
   * \param[in] val_consvar_grad_j - Gradient of the conservative variable at point j.
   */
  void SetConsVarGradient(su2double **val_consvar_grad_i, su2double **val_consvar_grad_j);
  
  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad_0 - Gradient of the conservative variable at point 0.
   * \param[in] val_consvar_grad_1 - Gradient of the conservative variable at point 1.
   * \param[in] val_consvar_grad_2 - Gradient of the conservative variable at point 2.
   */
  void SetConsVarGradient(su2double **val_consvar_grad_0,
                          su2double **val_consvar_grad_1,
                          su2double **val_consvar_grad_2);
  
  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad_0 - Gradient of the conservative variable at point 0.
   * \param[in] val_consvar_grad_1 - Gradient of the conservative variable at point 1.
   * \param[in] val_consvar_grad_2 - Gradient of the conservative variable at point 2.
   * \param[in] val_consvar_grad_3 - Gradient of the conservative variable at point 3.
   */
  void SetConsVarGradient(su2double **val_consvar_grad_0,
                          su2double **val_consvar_grad_1,
                          su2double **val_consvar_grad_2,
                          su2double **val_consvar_grad_3);
  
  /*!
   * \brief Set the gradient of the conservative variables.
   * \param[in] val_consvar_grad - Gradient of the conservative variable which is a scalar.
   */
  void SetConsVarGradient(su2double **val_consvar_grad);
  
  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
   * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
   */
  void SetPrimVarGradient(su2double **val_primvar_grad_i,
                          su2double **val_primvar_grad_j);
  
  /*!
   * \brief Set the Limiter of the primitive variables.
   * \param[in] val_primvar_lim_i - Limiter of the primitive variable at point i.
   * \param[in] val_primvar_lim_j - Limiter of the primitive variable at point j.
   */

  void SetPrimVarLimiter(su2double *val_primvar_lim_i,
                         su2double *val_primvar_lim_j);
  
  /*!
   * \brief Set the value of the adjoint variable.
   * \param[in] val_psi_i - Value of the adjoint variable at point i.
   * \param[in] val_psi_j - Value of the adjoint variable at point j.
   */
  void SetAdjointVar(su2double *val_psi_i, su2double *val_psi_j);
  
  /*!
   * \brief Set the value of the linearized conservative variables.
   * \param[in] val_deltau_i - Value of the linearized conservative variable at point i.
   * \param[in] val_deltau_j - Value of the linearized conservative variable at point j.
   */
  void SetLinearizedVar(su2double *val_deltau_i, su2double *val_deltau_j);
  
  /*!
   * \brief Set the gradient of the adjoint variables.
   * \param[in] val_psivar_grad_i - Gradient of the adjoint variable at point i.
   * \param[in] val_psivar_grad_j - Gradient of the adjoint variable at point j.
   */
  void SetAdjointVarGradient(su2double **val_psivar_grad_i, su2double **val_psivar_grad_j);
  
  /*!
   * \brief Set the limiter of the adjoint variables.
   * \param[in] val_psivar_lim_i - Gradient of the adjoint variable at point i.
   * \param[in] val_psivar_lim_j - Gradient of the adjoint variable at point j.
   */
  void SetAdjointVarLimiter(su2double *val_psivar_lim_i, su2double *val_psivar_lim_j);
  
  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
   * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
   */
  void SetTurbVar(su2double *val_turbvar_i, su2double *val_turbvar_j);
  
  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_transvar_i - Value of the turbulent variable at point i.
   * \param[in] val_transvar_j - Value of the turbulent variable at point j.
   */
  void SetTransVar(su2double *val_transvar_i, su2double *val_transvar_j);
  
  /*!
   * \brief Set the gradient of the turbulent variables.
   * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
   */
  void SetTurbVarGradient(su2double **val_turbvar_grad_i, su2double **val_turbvar_grad_j);
  
  /*!
   * \brief Set the gradient of the turbulent variables.
   * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
   */
  void SetTransVarGradient(su2double **val_transvar_grad_i, su2double **val_transvar_grad_j);
  
  /*!
   * \brief Set the value of the level set variable.
   * \param[in] val_levelsetvar_i - Value of the level set variable at point i.
   * \param[in] val_levelsetvar_j - Value of the level set variable at point j.
   */
  void SetLevelSetVar(su2double *val_levelsetvar_i, su2double *val_levelsetvar_j);
  
  /*!
   * \brief Set the gradient of the level set variables.
   * \param[in] val_levelsetvar_grad_i - Gradient of the level set variable at point i.
   * \param[in] val_levelsetvar_grad_j - Gradient of the level set variable at point j.
   */
  void SetLevelSetVarGradient(su2double **val_levelsetvar_grad_i, su2double **val_levelsetvar_grad_j);
  
  /*!
   * \brief Set the value of the adjoint turbulent variable.
   * \param[in] val_turbpsivar_i - Value of the adjoint turbulent variable at point i.
   * \param[in] val_turbpsivar_j - Value of the adjoint turbulent variable at point j.
   */
  void SetTurbAdjointVar(su2double *val_turbpsivar_i, su2double *val_turbpsivar_j);
  
  /*!
   * \brief Set the gradient of the adjoint turbulent variables.
   * \param[in] val_turbpsivar_grad_i - Gradient of the adjoint turbulent variable at point i.
   * \param[in] val_turbpsivar_grad_j - Gradient of the adjoint turbulent variable at point j.
   */
  void SetTurbAdjointGradient (su2double **val_turbpsivar_grad_i, su2double **val_turbpsivar_grad_j);
  
  /*!
   * \brief Set the value of the first blending function.
   * \param[in] val_F1_i - Value of the first Menter blending function at point i.
   * \param[in] val_F1_j - Value of the first Menter blending function at point j.
   */
  virtual void SetF1blending(su2double val_F1_i, su2double val_F1_j) {/* empty */};
  
  /*!
   * \brief Set the value of the second blending function.
   * \param[in] val_F1_i - Value of the second Menter blending function at point i.
   * \param[in] val_F1_j - Value of the second Menter blending function at point j.
   */
  virtual void SetF2blending(su2double val_F1_i, su2double val_F1_j) {/* empty */};
  
  /*!
   * \brief Set the value of the cross diffusion for the SST model.
   * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
   * \param[in] val_CDkw_j - Value of the cross diffusion at point j.
   */
  virtual void SetCrossDiff(su2double val_CDkw_i, su2double val_CDkw_j) {/* empty */};
  
  /*!
   * \brief Set the gradient of the auxiliary variables.
   * \param[in] val_auxvargrad_i - Gradient of the auxiliary variable at point i.
   * \param[in] val_auxvargrad_j - Gradient of the auxiliary variable at point j.
   */
  void SetAuxVarGrad(su2double *val_auxvargrad_i, su2double *val_auxvargrad_j);
  
  /*!
   * \brief Set the diffusion coefficient
   * \param[in] val_diffusioncoeff_i - Value of the diffusion coefficients at i.
   * \param[in] val_diffusioncoeff_j - Value of the diffusion coefficients at j
   */
  void SetDiffusionCoeff(su2double* val_diffusioncoeff_i,
                         su2double* val_diffusioncoeff_j);
  
  /*!
   * \brief Set the laminar viscosity.
   * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
   * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
   */
  void SetLaminarViscosity(su2double val_laminar_viscosity_i,
                           su2double val_laminar_viscosity_j);
  
  /*!
   * \brief Set the thermal conductivity (translational/rotational)
   * \param[in] val_thermal_conductivity_i - Value of the thermal conductivity at point i.
   * \param[in] val_thermal_conductivity_j - Value of the thermal conductivity at point j.
   * \param[in] iSpecies - Value of the species.
   */
  void SetThermalConductivity(su2double val_thermal_conductivity_i,
                              su2double val_thermal_conductivity_j);
  
  /*!
   * \brief Set the thermal conductivity (translational/rotational)
   * \param[in] val_thermal_conductivity_i - Value of the thermal conductivity at point i.
   * \param[in] val_thermal_conductivity_j - Value of the thermal conductivity at point j.
   * \param[in] iSpecies - Value of the species.
   */
  void SetThermalConductivity_ve(su2double val_thermal_conductivity_ve_i,
                                 su2double val_thermal_conductivity_ve_j);
  
  /*!
   * \brief Set the eddy viscosity.
   * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
   * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
   */
  void SetEddyViscosity(su2double val_eddy_viscosity_i,
                        su2double val_eddy_viscosity_j);
  
  /*!
   * \brief Set the turbulent kinetic energy.
   * \param[in] val_turb_ke_i - Value of the turbulent kinetic energy at point i.
   * \param[in] val_turb_ke_j - Value of the turbulent kinetic energy at point j.
   */
  void SetTurbKineticEnergy(su2double val_turb_ke_i, su2double val_turb_ke_j);
  
  /*!
   * \brief Set the value of the distance from the nearest wall.
   * \param[in] val_dist_i - Value of of the distance from point i to the nearest wall.
   * \param[in] val_dist_j - Value of of the distance from point j to the nearest wall.
   */
  void SetDistance(su2double val_dist_i, su2double val_dist_j);
  
  /*!
   * \brief Set coordinates of the points.
   * \param[in] val_coord_i - Coordinates of the point i.
   * \param[in] val_coord_j - Coordinates of the point j.
   */
  void SetCoord(su2double *val_coord_i, su2double *val_coord_j);
  
  /*!
   * \overload
   * \param[in] val_coord_0 - Coordinates of the point 0.
   * \param[in] val_coord_1 - Coordinates of the point 1.
   * \param[in] val_coord_2 - Coordinates of the point 2.
   */
  void SetCoord(su2double *val_coord_0, su2double *val_coord_1, su2double *val_coord_2);
  
  /*!
   * \overload
   * \param[in] val_coord_0 - Coordinates of the point 0.
   * \param[in] val_coord_1 - Coordinates of the point 1.
   * \param[in] val_coord_2 - Coordinates of the point 2.
   * \param[in] val_coord_3 - Coordinates of the point 3.
   */
  void SetCoord(su2double *val_coord_0, su2double *val_coord_1, su2double *val_coord_2,
                su2double *val_coord_3);
  
  /*!
   * \brief Set the velocity of the computational grid.
   * \param[in] val_gridvel_i - Grid velocity of the point i.
   * \param[in] val_gridvel_j - Grid velocity of the point j.
   */
  void SetGridVel(su2double *val_gridvel_i, su2double *val_gridvel_j);
  
  /*!
   * \brief Set the wind gust value.
   * \param[in] val_windgust_i - Wind gust of the point i.
   * \param[in] val_windgust_j - Wind gust of the point j.
   */
  void SetWindGust(su2double *val_windgust_i, su2double *val_windgust_j);
  
  /*!
   * \brief Set the wind gust derivatives values.
   * \param[in] val_windgust_i - Wind gust derivatives of the point i.
   * \param[in] val_windgust_j - Wind gust derivatives of the point j.
   */
  void SetWindGustDer(su2double *val_windgustder_i, su2double *val_windgustder_j);
  
  /*!
   * \brief Set the value of the pressure.
   * \param[in] val_pressure_i - Value of the pressure at point i.
   * \param[in] val_pressure_j - Value of the pressure at point j.
   */
  void SetPressure(su2double val_pressure_i, su2double val_pressure_j);
  
  /*!
   * \brief Set the value of the density for the incompressible solver.
   * \param[in] val_densityinc_i - Value of the pressure at point i.
   * \param[in] val_densityinc_j - Value of the pressure at point j.
   */
  void SetDensityInc(su2double val_densityinc_i, su2double val_densityinc_j);
  
  /*!
   * \brief Set the value of the beta for incompressible flows.
   * \param[in] val_betainc2_i - Value of beta for incompressible flows at point i.
   * \param[in] val_betainc2_j - Value of beta for incompressible flows at point j.
   */
  void SetBetaInc2(su2double val_betainc2_i, su2double val_betainc2_j);
  
  /*!
   * \brief Set the value of the sound speed.
   * \param[in] val_soundspeed_i - Value of the sound speed at point i.
   * \param[in] val_soundspeed_j - Value of the sound speed at point j.
   */
  void SetSoundSpeed(su2double val_soundspeed_i, su2double val_soundspeed_j);
  
  /*!
   * \brief Set the value of the temperature.
   * \param[in] val_temp_i - Value of the temperature at point i.
   * \param[in] val_temp_j - Value of the temperature at point j.
   */
  void SetTemperature(su2double val_temp_i, su2double val_temp_j);
  
  /*!
   * \brief Set the value of the enthalpy.
   * \param[in] val_enthalpy_i - Value of the enthalpy at point i.
   * \param[in] val_enthalpy_j - Value of the enthalpy at point j.
   */
  void SetEnthalpy(su2double val_enthalpy_i, su2double val_enthalpy_j);
  
  /*!
   * \brief Set the value of the spectral radius.
   * \param[in] val_lambda_i - Value of the spectral radius at point i.
   * \param[in] val_lambda_j - Value of the spectral radius at point j.
   */
  void SetLambda(su2double val_lambda_i, su2double val_lambda_j);
  
  /*!
   * \brief Set the value of undivided laplacian.
   * \param[in] val_und_lapl_i Undivided laplacian at point i.
   * \param[in] val_und_lapl_j Undivided laplacian at point j.
   */
  void SetUndivided_Laplacian(su2double *val_und_lapl_i, su2double *val_und_lapl_j);
  
  /*!
   * \brief Set the value of the pressure sensor.
   * \param[in] val_sensor_i Pressure sensor at point i.
   * \param[in] val_sensor_j Pressure sensor at point j.
   */
  void SetSensor(su2double val_sensor_i, su2double val_sensor_j);
  
  /*!
   * \brief Set the number of neighbor to a point.
   * \param[in] val_neighbor_i - Number of neighbor to point i.
   * \param[in] val_neighbor_j - Number of neighbor to point j.
   */
  void SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j);
  
  /*!
   * \brief Set the value of the normal vector to the face between two points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void SetNormal(su2double *val_normal);
  
  /*!
   * \brief Set the value of the volume of the control volume.
   * \param[in] val_volume Volume of the control volume.
   */
  void SetVolume(su2double val_volume);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetRhosIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetRhoIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetPIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetTIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetTveIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the velocity index in the primitive variable vector.
   * \param[in] i(rho*u)
   */
  void SetVelIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetHIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetAIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetRhoCvtrIndex(unsigned short val_Index);
  
  /*!
   * \brief Retrieves the value of the species density in the primitive variable vector.
   * \param[in] iRho_s
   */
  void SetRhoCvveIndex(unsigned short val_Index);
  
  /*!
   * \brief Sets the value of the derivative of pressure w.r.t. species density.
   * \param[in] iRho_s
   */
  void SetdPdU(su2double *val_dPdU_i, su2double *val_dPdU_j);
  
  /*!
   * \brief Sets the value of the derivative of temperature w.r.t. species density.
   * \param[in] iRho_s
   */
  void SetdTdU(su2double *val_dTdU_i, su2double *val_dTdU_j);
  
  /*!
   * \brief Sets the value of the derivative of vib-el. temperature w.r.t. species density.
   * \param[in] iRho_s
   */
  void SetdTvedU(su2double *val_dTvedU_i, su2double *val_dTvedU_j);
  
  /*!
   * \brief Get the inviscid fluxes.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_pressure - Value of the pressure.
   * \param[in] val_enthalpy - Value of the enthalpy.
   */
  void GetInviscidFlux(su2double val_density, su2double *val_velocity, su2double val_pressure, su2double val_enthalpy);
  
  /*!
   * \brief Get the viscous fluxes.
   * \param[in] val_primvar - Value of the primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_mach_inf - Value of the Mach number at the infinity.
   */
  void GetViscousFlux(su2double *val_primvar, su2double **val_gradprimvar,
                      su2double val_laminar_viscosity, su2double val_eddy_viscosity,
                      su2double val_mach_inf);
  
  /*!
   * \brief Compute the projected inviscid flux vector.
   * \param[in] val_density - Pointer to the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_pressure - Pointer to the pressure.
   * \param[in] val_enthalpy - Pointer to the enthalpy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidProjFlux(su2double *val_density, su2double *val_velocity,
                           su2double *val_pressure, su2double *val_enthalpy,
                           su2double *val_normal, su2double *val_Proj_Flux);
  
  /*!
   * \brief Compute the projected inviscid flux vector.
   * \param[in] val_U - Conserved variables
   * \param[in] val_V - Primitive variables
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidProjFlux(su2double *val_U, su2double *val_V, su2double *val_normal,
                           su2double *val_Proj_Flux);
  
  /*!
   * \brief Compute the projected inviscid flux vector for incompresible simulations
   * \param[in] val_density - Pointer to the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_pressure - Pointer to the pressure.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidArtCompProjFlux(su2double *val_density, su2double *val_velocity,
                                  su2double *val_pressure, su2double *val_betainc2,
                                  su2double *val_normal, su2double *val_Proj_Flux);
  
  /*!
   * \brief Compute the projected inviscid flux vector for incompresible simulations
   * \param[in] val_density - Pointer to the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_pressure - Pointer to the pressure.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_Proj_Flux - Pointer to the projected flux.
   */
  void GetInviscidArtComp_FreeSurf_ProjFlux(su2double *val_density,
                                            su2double *val_velocity,
                                            su2double *val_pressure,
                                            su2double *val_betainc2,
                                            su2double *val_levelset,
                                            su2double *val_normal,
                                            su2double *val_Proj_Flux);
  
  
  /*!
   * \brief Compute the projection of the viscous fluxes into a direction.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_turb_ke - Turbulent kinetic energy
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_eddy_conductivity - Eddy Conductivity.
   */
  
  void GetViscousProjFlux(su2double *val_primvar, su2double **val_gradprimvar,
                          su2double val_turb_ke, su2double *val_normal,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity);
  /*!
   * \brief Compute the projection of the viscous fluxes into a direction for general fluid model.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_turb_ke - Turbulent kinetic energy
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_heat_capacity_cp - Heat Capacity at constant pressure.
   */
  
  void GetViscousProjFlux(su2double *val_primvar, su2double **val_gradprimvar,
                          su2double val_turb_ke, su2double *val_normal,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity,
                          su2double val_thermal_conductivity,
                          su2double val_heat_capacity_cp);
  
  
  
  /*!
   * * \brief Compute the projection of the viscous fluxes into a direction.
   * \brief Overloaded function for multiple species viscous calculations
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_diffusioncoeff
   * \param[in] val_therm_conductivity
   * \param[in] val_therm_conductivity_ve
   * \param[in] config
   */
  void GetViscousProjFlux(su2double *val_primvar,
                          su2double **val_gradprimvar,
                          su2double *val_normal,
                          su2double *val_diffusioncoeff,
                          su2double val_viscosity,
                          su2double val_therm_conductivity,
                          su2double val_therm_conductivity_ve,
                          CConfig *config);
  
  /*
   * \brief Compute the projection of the viscous fluxes into a direction (artificial compresibility method).
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   */
  
  void GetViscousArtCompProjFlux(su2double **val_gradprimvar,
                                 su2double *val_normal,
                                 su2double val_laminar_viscosity,
                                 su2double val_eddy_viscosity);
  
  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices.
   * \param[in] val_velocity Pointer to the velocity.
   * \param[in] val_energy Value of the energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(su2double *val_velocity, su2double *val_energy,
                          su2double *val_normal, su2double val_scale,
                          su2double **val_Proj_Jac_tensor);
  
  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices (artificial compresibility).
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidArtCompProjJac(su2double *val_density, su2double *val_velocity,
                                 su2double *val_betainc2, su2double *val_normal,
                                 su2double val_scale,
                                 su2double **val_Proj_Jac_tensor);
  
  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices (artificial compresibility).
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Pointer to the velocity.
   * \param[in] val_betainc2 - Value of the artificial compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidArtComp_FreeSurf_ProjJac(su2double *val_density,
                                           su2double *val_ddensity,
                                           su2double *val_velocity,
                                           su2double *val_betainc2,
                                           su2double *val_levelset,
                                           su2double *val_normal,
                                           su2double val_scale,
                                           su2double **val_Proj_Jac_tensor);
  
  /*!
   * \brief Compute the projection of the inviscid Jacobian matrices for general fluid model.
   * \param[in] val_velocity Pointer to the velocity.
   * \param[in] val_energy Value of the energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(su2double *val_velocity, su2double *val_enthalphy,
                          su2double *val_chi, su2double *val_kappa,
                          su2double *val_normal, su2double val_scale,
                          su2double **val_Proj_Jac_tensor);
  
  /*!
   * \overload
   * \brief Compute the projection of the inviscid Jacobian matrices for the two-temperature model.
   * \param[in] val_U - Vector conserved variables.
   * \param[in] val_V - Vector of primitive variables.
   * \param[in] val_dPdU - Vector of partial derivatives of pressure w.r.t. conserved vars.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_scale - Scale of the projection.
   * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
   */
  void GetInviscidProjJac(su2double *val_U, su2double *val_V, su2double *val_dPdU,
                          su2double *val_normal, su2double val_scale,
                          su2double **val_Proj_Jac_Tensor);
  
  /*!
   * \brief TSL-Approximation of Viscous NS Jacobians.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(su2double *val_Mean_PrimVar,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity,
                          su2double val_dist_ij,
                          su2double *val_normal, su2double val_dS,
                          su2double *val_Proj_Visc_Flux,
                          su2double **val_Proj_Jac_Tensor_i,
                          su2double **val_Proj_Jac_Tensor_j);
  
  /*!
   * \brief TSL-Approximation of Viscous NS Jacobians for arbitrary equations of state.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_gradprimvar - Mean value of the gradient of the primitive variables.
   * \param[in] val_Mean_SecVar - Mean value of the secondary variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_heat_capacity_cp - Value of the specific heat at constant pressure.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(su2double *val_Mean_PrimVar,
                          su2double **val_gradprimvar,
                          su2double *val_Mean_SecVar,
                          su2double val_laminar_viscosity,
                          su2double val_eddy_viscosity,
                          su2double val_thermal_conductivity,
                          su2double val_heat_capacity_cp,
                          su2double val_dist_ij,
                          su2double *val_normal, su2double val_dS,
                          su2double *val_Proj_Visc_Flux,
                          su2double **val_Proj_Jac_Tensor_i,
                          su2double **val_Proj_Jac_Tensor_j);
  
  /*!
   * \brief Mapping between primitives variables P and conservatives variables C.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_Mean_PrimVar - Mean Value of the secondary variables.
   * \param[out] val_Jac_PC - Pointer to the Jacobian dPdC.
   */
  void GetPrimitive2Conservative (su2double *val_Mean_PrimVar,
                                  su2double *val_Mean_SecVar,
                                  su2double **val_Jac_PC);
  
  /*!
   * \brief TSL-Approximation of Viscous NS Jacobians.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_thermal_conductivity
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(su2double *val_Mean_PrimVar,
                          su2double *val_diffusion_coeff,
                          su2double val_laminar_viscosity,
                          su2double val_thermal_conductivity,
                          su2double val_thermal_conductivity_ve,
                          su2double val_dist_ij,
                          su2double *val_normal, su2double val_dS,
                          su2double *val_Proj_Visc_Flux,
                          su2double **val_Proj_Jac_Tensor_i,
                          su2double **val_Proj_Jac_Tensor_j,
                          CConfig *config);
  
  /*!
   * \brief Compute the projection of the viscous Jacobian matrices.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousArtCompProjJacs(su2double val_laminar_viscosity,
                                 su2double val_eddy_viscosity, su2double val_dist_ij,
                                 su2double *val_normal, su2double val_dS,
                                 su2double **val_Proj_Jac_Tensor_i,
                                 su2double **val_Proj_Jac_Tensor_j);
  
  /*!
   * \overload
   * \brief Computation of the matrix P for a generic fluid model
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_enthalpy - Value of the Enthalpy
   * \param[in] val_chi - Value of the derivative of Pressure with respect to the Density.
   * \param[in] val_kappa - Value of the derivative of Pressure with respect to the volume specific Static Energy.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_p_tensor - Pointer to the P matrix.
   */
  void GetPMatrix(su2double *val_density, su2double *val_velocity,
                  su2double *val_soundspeed, su2double *val_enthalpy, su2double *val_chi, su2double *val_kappa,
                  su2double *val_normal, su2double **val_p_tensor);
  
  /*!
   * \brief Computation of the matrix P, this matrix diagonalize the conservative Jacobians in
   *        the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_p_tensor - Pointer to the P matrix.
   */
  void GetPMatrix(su2double *val_density, su2double *val_velocity,
                  su2double *val_soundspeed, su2double *val_normal,
                  su2double **val_p_tensor);
  
  /*!
   * \overload
   * \brief Computation of the matrix P, this matrix diagonalizes the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] U - Vector of conserved variables (really only need rhoEve)
   * \param[in] V - Vector of primitive variables
   * \param[in] val_dPdU - Vector of derivatives of pressure w.r.t. conserved vars.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] l - Tangential vector to face.
   * \param[in] m - Tangential vector to face (mutually orthogonal to val_normal & l).
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix(su2double *U, su2double *V, su2double *val_dPdU,
                  su2double *val_normal, su2double *l, su2double *m,
                  su2double **val_p_tensor) ;
  
  /*!
   * \brief Computation of the matrix Rinv*Pe.
   * \param[in] Beta2 - A variable in used to define Pe matrix.
   * \param[in] val_enthalpy - value of the enthalpy.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[in] val_velocity - value of the velocity.
   * \param[out] val_invR_invPe - Pointer to the matrix of conversion from entropic to conserved variables.
   */
  void GetinvRinvPe(su2double Beta2, su2double val_enthalpy, su2double val_soundspeed,
                    su2double val_density, su2double* val_velocity,
                    su2double** val_invR_invPe);
  
  /*!
   * \brief Computation of the matrix R.
   * \param[in] val_pressure - value of the pressure.
   * \param[in] val_soundspeed - value of the sound speed.
   * \param[in] val_density - value of the density.
   * \param[in] val_velocity - value of the velocity.
   * \param[out] val_invR_invPe - Pointer to the matrix of conversion from entropic to conserved variables.
   */
  void GetRMatrix(su2double val_pressure, su2double val_soundspeed,
                  su2double val_density, su2double* val_velocity,
                  su2double** val_invR_invPe);
  
  /*!
   * \brief Computation of the matrix Td, this matrix diagonalize the preconditioned conservative Jacobians
   *        in the form $Tg |Lambda| Td = Pc{-1}|Pc (A.Normal)|$.
   * \param[in] Beta2 - A variable in used to define absPeJacobian matrix.
   * \param[in] r_hat - A variable in used to define absPeJacobian matrix.
   * \param[in] s_hat - A variable in used to define absPeJacobian matrix.
   * \param[in] t_hat - A variable in used to define absPeJacobian matrix.
   * \param[in] rB2a2 - A variable in used to define absPeJacobian matrix.
   * \param[in] val_Lambda - Eigenvalues of the Preconditioned Jacobian.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_absPeJac - Pointer to the Preconditioned Jacobian matrix.
   */
  void GetPrecondJacobian(su2double Beta2, su2double r_hat, su2double s_hat, su2double t_hat, su2double rB2a2, su2double* val_Lambda, su2double* val_normal, su2double** val_absPeJac);
  
  /*!
   * \brief Computation of the matrix P (artificial compresibility), this matrix diagonalize the conservative Jacobians in
   *        the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_betainv2 - Value of the compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_p_tensor - Pointer to the P matrix.
   */
  void GetPArtCompMatrix(su2double *val_density, su2double *val_velocity,
                         su2double *val_betainv2, su2double *val_normal,
                         su2double **val_p_tensor);
  
  /*!
   * \brief Computation of the matrix P (artificial compresibility), this matrix diagonalize the conservative Jacobians in
   *        the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_betainv2 - Value of the compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_p_tensor - Pointer to the P matrix.
   */
  void GetPArtComp_FreeSurf_Matrix(su2double *val_density, su2double *val_ddensity,
                                   su2double *val_velocity, su2double *val_betainv2,
                                   su2double *val_levelset, su2double *val_normal,
                                   su2double **val_p_tensor);
  
  /*!
   * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians
   * in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(su2double **val_invp_tensor, su2double *val_density,
                      su2double *val_velocity, su2double *val_soundspeed,
                      su2double *val_chi, su2double *val_kappa,
                      su2double *val_normal);
  
  /*!
   * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(su2double *val_density, su2double *val_velocity,
                      su2double *val_soundspeed, su2double *val_normal,
                      su2double **val_invp_tensor);
  
  /*!
   * \overload
   * \brief Computation of the matrix P^{-1}, this matrix diagonalizes the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] U - Vector of conserved variables.
   * \param[in] V - Vector of primitive variables.
   * \param[in] val_dPdU - Vector of derivatives of pressure w.r.t. conserved variables
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] l - Tangential vector to face.
   * \param[in] m - Tangential vector to face (mutually orthogonal to val_normal & l).
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPMatrix_inv(su2double *U, su2double *V, su2double *val_dPdU,
                      su2double *val_normal, su2double *l, su2double *m,
                      su2double **val_invp_tensor) ;
  
  /*!
   * \brief Computation of the matrix P^{-1} (artificial compresibility), this matrix diagonalize the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_betainv2 - Value of the compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPArtCompMatrix_inv(su2double *val_density, su2double *val_velocity,
                             su2double *val_betainv2, su2double *val_normal,
                             su2double **val_invp_tensor);
  
  /*!
   * \brief Computation of the matrix P^{-1} (artificial compresibility), this matrix diagonalize the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] val_density - Value of the density.
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_betainv2 - Value of the compresibility factor.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
   */
  void GetPArtComp_FreeSurf_Matrix_inv(su2double *val_density,
                                       su2double *val_ddensity,
                                       su2double *val_velocity,
                                       su2double *val_betainv2,
                                       su2double *val_levelset,
                                       su2double *val_normal,
                                       su2double **val_invp_tensor);
  
  /*!
   * \brief Compute viscous residual and jacobian.
   */
  void GetAdjViscousFlux_Jac(su2double Pressure_i, su2double Pressure_j, su2double Density_i, su2double Density_j,
                             su2double ViscDens_i, su2double ViscDens_j, su2double *Velocity_i, su2double *Velocity_j,
                             su2double sq_vel_i, su2double sq_vel_j,
                             su2double XiDens_i, su2double XiDens_j, su2double **Mean_GradPhi, su2double *Mean_GradPsiE,
                             su2double dPhiE_dn, su2double *Normal, su2double *Edge_Vector, su2double dist_ij_2, su2double *val_residual_i,
                             su2double *val_residual_j,
                             su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji,
                             su2double **val_Jacobian_jj, bool implicit);
  
  /*!
   * \brief Computation of the projected inviscid lambda (eingenvalues).
   * \param[in] val_velocity - Value of the velocity.
   * \param[in] val_soundspeed - Value of the sound speed.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_Lambda_Vector - Pointer to Lambda matrix.
   */
  void GetJacInviscidLambda_fabs(su2double *val_velocity, su2double val_soundspeed,
                                 su2double *val_normal, su2double *val_Lambda_Vector);
  
  /*!
   * \brief Compute the numerical residual.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_residual, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  virtual void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j);
  
  virtual void ComputeResidual_TransLM(su2double *val_residual,
                                       su2double **val_Jacobian_i,
                                       su2double **val_Jacobian_j, CConfig *config,
                                       su2double &gamma_sep) ;
  
  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_residual_i,
                               su2double *val_residual_j, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                               su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[out] val_JacobianMeanFlow_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_JacobianMeanFlow_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                               su2double **val_Jacobian_j,
                               su2double **val_JacobianMeanFlow_i,
                               su2double **val_JacobianMeanFlow_j,
                               CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                               CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_resconv, su2double *val_resvisc,
                               su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                               CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                               su2double **val_Jacobian_ii,
                               su2double **val_Jacobian_ij,
                               su2double **val_Jacobian_ji,
                               su2double **val_Jacobian_jj, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double *val_resconv_i, su2double *val_resvisc_i,
                               su2double *val_resconv_j, su2double *val_resvisc_j,
                               su2double **val_Jacobian_ii,
                               su2double **val_Jacobian_ij,
                               su2double **val_Jacobian_ji,
                               su2double **val_Jacobian_jj, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double **val_stiffmatrix_elem, CConfig *config);
  
  /*!
   * \overload
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  virtual void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                               CConfig *config);
  
  /*!
   * \overload
   * \param[out] - Matrix for storing the constants to be used in the calculation of the equilibrium extent of reaction Keq.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GetEq_Rxn_Coefficients(su2double **EqnRxnConstants, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the source residual containing chemistry terms.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual_Axisymmetric(su2double *val_residual, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the source residual containing chemistry terms.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual_Axisymmetric_ad(su2double *val_residual, su2double *val_residuald, CConfig *config);
  
  /*!
   * \brief Calculation of axisymmetric source term Jacobian
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetJacobian_Axisymmetric(su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Calculation of the translational-vibrational energy exchange source term
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  virtual void ComputeVibRelaxation(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Calculation of the chemistry source term
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   * \param[out] val_Jacobian_i - Jacobian of the source terms
   */
  virtual void ComputeChemistry(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Calculates constants used for Keq correlation.
   * \param[out] A - Pointer to coefficient array.
   * \param[in] val_reaction - Reaction number indicator.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GetKeqConstants(su2double *A, unsigned short val_reaction, CConfig *config);
  
  /*!
   * \brief Set intermittency for numerics (used in SA with LM transition model)
   */
  virtual void SetIntermittency(su2double intermittency_in);
  
  /*!
   * \brief Computes the viscous source term for the TNE2 adjoint problem
   * \param[in] config - Definition of the particular problem.
   * \param[out] val_residual - residual of the source terms
   */
  virtual void ComputeSourceViscous(su2double *val_residual, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  virtual void SetProduction(su2double val_production);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  virtual void SetDestruction(su2double val_destruction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  virtual void SetCrossProduction(su2double val_crossproduction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  virtual su2double GetProduction(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  virtual su2double GetDestruction(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  virtual su2double GetCrossProduction(void);
  
  /*!
   * \overload
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double **val_Jacobian_i,
                               su2double *val_Jacobian_mui,
                               su2double ***val_Jacobian_gradi, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeResidual(su2double **val_Jacobian_i,
                               su2double *val_Jacobian_mui,
                               su2double ***val_Jacobian_gradi,
                               su2double **val_Jacobian_j,
                               su2double *val_Jacobian_muj,
                               su2double ***val_Jacobian_gradj, CConfig *config);
  
  /*!
   * \brief Computing stiffness matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetFEA_StiffMatrix2D(su2double **StiffMatrix_Elem, su2double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d);
  
  /*!
   * \brief Computing stiffness matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetFEA_StiffMatrix3D(su2double **StiffMatrix_Elem, su2double CoordCorners[8][3], unsigned short nNodes);
  
  /*!
   * \brief Computing mass matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetFEA_StiffMassMatrix2D(su2double **StiffMatrix_Elem, su2double **MassMatrix_Elem, su2double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d);
  
  /*!
   * \brief Computing mass matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetFEA_StiffMassMatrix3D(su2double **StiffMatrix_Elem, su2double **MassMatrix_Elem, su2double CoordCorners[8][3], unsigned short nNodes);
  
  /*!
   * \brief Computing dead load vector of the Galerkin method.
   * \param[out] val_deadloadvector_elem - Dead load at the nodes for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetFEA_DeadLoad2D(su2double *DeadLoadVector_Elem, su2double CoordCorners[8][3], unsigned short nNodes, su2double matDensity);
  
  /*!
   * \brief Computing stiffness matrix of the Galerkin method.
   * \param[out] val_deadloadvector_elem - Dead load at the nodes for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetFEA_DeadLoad3D(su2double *DeadLoadVector_Elem, su2double CoordCorners[8][3], unsigned short nNodes, su2double matDensity);
  
  
  /*!
   * \brief Computing stresses in FEA method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GetFEA_StressNodal2D(su2double StressVector[8][3], su2double DispElement[8], su2double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d);
  
  
  /*!
   * \brief Computing stresses in FEA method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void GetFEA_StressNodal3D(su2double StressVector[8][6], su2double DispElement[24], su2double CoordCorners[8][3], unsigned short nNodes);
  
  /*!
   * \brief A virtual member to linearly interpolate pressures
   * \param[in] config - Definition of the particular problem.
   */
  virtual void PressInt_Linear(su2double CoordCorners[4][3], su2double *tn_e, su2double Fnodal[12]);
  
  /*!
   * \brief A virtual member to linearly interpolate viscous stresses
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ViscTermInt_Linear(su2double CoordCorners[2][2], su2double Tau_0[3][3], su2double Tau_1[3][3],  su2double FviscNodal[4]);
  
  /*!
   * \brief Computes a basis of orthogonal vectors from a suppled vector
   * \param[in] config - Normal vector
   */
  void CreateBasis(su2double *val_Normal);
  
};

/*!
 * \class CUpwCUSP_Flow
 * \brief Class for centered scheme - CUSP.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CUpwCUSP_Flow : public CNumerics {
  
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, *Diff_Flux, /*!< \brief Diference of conservative variables and undivided laplacians. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity,  /*!< \brief Mean and projected velocities. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  *ProjFlux, *ProjFlux_i, *ProjFlux_j,  /*!< \brief Projected inviscid flux tensor. */
  cte_0, cte_1, /*!< \brief Artificial dissipation values. */
  LamdaNeg, LamdaPos, Beta, Nu_c, U_i[5], U_j[5], MeanSoundSpeed, Mach,
  **Jacobian;  /*!< \brief Projected grid velocity. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwCUSP_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwCUSP_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_residual - Pointer to the residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CUpwRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations.
 * \ingroup ConvDiscr
 * \author A. Bueno, F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CUpwRoe_Flow : public CNumerics {
private:
  bool implicit, grid_movement;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon, MaxLambda, Delta;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho, RoeSoundSpeed2, kappa;
  unsigned short iDim, iVar, jVar, kVar;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};


/*!
 * \class CUpwGeneralRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations for a general fluid model.
 * \ingroup ConvDiscr
 * \author S.Vitale, G.Gori, M.Pini
 * \version 4.0.1 "Cardinal"
 */
class CUpwGeneralRoe_Flow : public CNumerics {
private:
  bool implicit, grid_movement;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon, MaxLambda, Delta;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed, RoeSoundSpeed2,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho, kappa;
  unsigned short iDim, iVar, jVar, kVar;
  
  
  su2double StaticEnthalpy_i, StaticEnergy_i, StaticEnthalpy_j, StaticEnergy_j, Kappa_i, Kappa_j, Chi_i, Chi_j, Velocity2_i, Velocity2_j;
  su2double RoeKappa, RoeChi, RoeKappaStaticEnthalpy;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwGeneralRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwGeneralRoe_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Compute the Average for a general fluid flux between two nodes i and j.
   * Using the approach of Vinokur and Montagne'
   */
  
  void ComputeRoeAverage();
};


/*!
 * \class CUpwMSW_Flow
 * \brief Class for solving a flux-vector splitting method by Steger & Warming, modified version.
 * \ingroup ConvDiscr
 * \author S. Copeland
 * \version 4.0.1 "Cardinal"
 */
class CUpwMSW_Flow : public CNumerics {
private:
  bool implicit;
  su2double *Diff_U;
  su2double *u_i, *u_j, *ust_i, *ust_j;
  su2double *Fc_i, *Fc_j;
  su2double *Lambda_i, *Lambda_j;
  su2double rhos_i, rhos_j;
  su2double *Ust_i, *Ust_j, *Vst_i, *Vst_j, *Velst_i, *Velst_j;
  su2double **P_Tensor, **invP_Tensor;
  unsigned short nPrimVar, nVar, nDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwMSW_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwMSW_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
};

/*!
 * \class CUpwTurkel_Flow
 * \brief Class for solving an approximate Riemann solver of Roe with Turkel Preconditioning for the flow equations.
 * \ingroup ConvDiscr
 * \author A. K. Lonkar
 * \version 4.0.1 "Cardinal"
 */
class CUpwTurkel_Flow : public CNumerics {
private:
  bool implicit, grid_movement;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **absPeJac, **invRinvPe, **R_Tensor, **Matrix, **Art_Visc;
  su2double sq_vel, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoePressure, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  unsigned short iDim, iVar, jVar, kVar;
  su2double Beta, Beta_min, Beta_max;
  su2double r_hat, s_hat, t_hat, rhoB2a2, sqr_one_m_Betasqr_Lam1;
  su2double Beta2, one_m_Betasqr, one_p_Betasqr, sqr_two_Beta_c_Area;
  su2double local_Mach;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwTurkel_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Get the Preconditioning Beta.
   * \return Beta - Value of the low Mach Preconditioner.
   */
  su2double GetPrecond_Beta();
};

/*!
 * \class CUpwArtComp_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the incompressible flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CUpwArtComp_Flow : public CNumerics {
private:
  bool implicit;
  bool gravity;
  su2double Froude;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *MeanVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double Proj_ModJac_Tensor_ij, Pressure_i,
  Pressure_j, MeanDensity, MeanSoundSpeed, MeanPressure, MeanBetaInc2,
  ProjVelocity;
  unsigned short iDim, iVar, jVar, kVar;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwArtComp_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwArtComp_FreeSurf_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the incompressible flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CUpwArtComp_FreeSurf_Flow : public CNumerics {
private:
  bool implicit;
  bool gravity;
  su2double Froude;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *MeanVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double Proj_ModJac_Tensor_ij, Pressure_i, LevelSet_i, dDensityInc_i, dDensityInc_j,
  Pressure_j, LevelSet_j, MeanDensityInc, dMeanDensityInc, MeanPressure, MeanLevelSet, MeanBetaInc2,
  ProjVelocity, Distance_i, Distance_j;
  unsigned short iDim, iVar, jVar, kVar;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwArtComp_FreeSurf_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwArtComp_FreeSurf_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwRoe_AdjFlow
 * \brief Class for solving an approximate Riemann solver of Roe
 *        for the adjoint flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CUpwRoe_AdjFlow : public CNumerics {
private:
  su2double *Residual_Roe;
  su2double area, Sx, Sy, Sz, rarea, nx, ny, nz, rho_l, u_l, v_l, w_l, h_l, rho_r,
  u_r, v_r, w_r, h_r, psi1, psi2, psi3, psi4, psi5;
  su2double h, u, v, w, c, psi1_l, psi2_l, psi3_l, psi4_l, psi5_l,
  psi1_r, psi2_r, psi3_r, psi4_r, psi5_r, q_l, q_r, Q_l, Q_r, vn,
  rrho_l, weight, rweight1, cc;
  su2double l1psi, l2psi, absQ, absQp, absQm, q2, alpha, beta_u, beta_v, beta_w, Q, l1l2p, l1l2m, eta;
  su2double RoeDensity, RoeSoundSpeed, *RoeVelocity, *Lambda, *Velocity_i, *Velocity_j, **ProjFlux_i, **ProjFlux_j,
  Proj_ModJac_Tensor_ij, **Proj_ModJac_Tensor, Energy_i, Energy_j, **P_Tensor, **invP_Tensor;
  unsigned short iDim, iVar, jVar, kVar;
  bool implicit, grid_movement;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoe_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint Roe's flux between two nodes i and j.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                       su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CUpwRoeArtComp_AdjFlow
 * \brief Class for solving an approximate Riemann solver of Roe
 *        for the adjoint flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CUpwRoeArtComp_AdjFlow : public CNumerics {
private:
  su2double Area, *Lambda, *Velocity_i, *Velocity_j, **Proj_Jac_Tensor_i, **Proj_Jac_Tensor_j,
  Proj_ModJac_Tensor_ij, **Proj_ModJac_Tensor, **P_Tensor, **invP_Tensor, MeanDensity,
  MeanPressure, MeanBetaInc2, ProjVelocity, *MeanVelocity, MeanSoundSpeed;
  unsigned short iDim, iVar, jVar, kVar;
  bool implicit;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoeArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoeArtComp_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint Roe's flux between two nodes i and j.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                       su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CUpwAUSM_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CUpwAUSM_Flow : public CNumerics {
private:
  bool implicit;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  unsigned short iDim, iVar, jVar, kVar;
  su2double mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSM_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwHLLC_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios, based on the Joe code implementation
 * \version 4.0.1 "Cardinal"
 */
class CUpwHLLC_Flow : public CNumerics {
private:
  bool implicit;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, sq_vel_i, sq_vel_j, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  unsigned short iDim, iVar, jVar, kVar;
  su2double Rrho, tmp, velRoe[3], uRoe, gamPdivRho, sq_velRoe, cRoe, sL, sR, sM, pStar, invSLmSs, sLmuL, rhoSL, rhouSL[3],
  eSL, invSRmSs, sRmuR, rhoSR, rhouSR[3], eSR;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwHLLC_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwLin_TransLM
 * \brief Class for performing a linear upwind solver for the Spalart-Allmaras turbulence model equations with transition
 * \ingroup ConvDiscr
 * \author A. Aranake
 * \version 4.0.1 "Cardinal"
 */
class CUpwLin_TransLM : public CNumerics {
private:
  su2double *Velocity_i;
  su2double *Velocity_j;
  bool implicit, grid_movement, incompressible;
  su2double Density_i, Density_j, q_ij, a0, a1;
  unsigned short iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwLin_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwLin_TransLM(void);
  
  /*!
   * \brief Compute the upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwLin_AdjLevelSet
 * \brief Class for performing a linear upwind solver for the adjoint Level Set equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CUpwLin_AdjLevelSet : public CNumerics {
private:
  bool implicit;
  su2double *Velocity_i;
  su2double *Velocity_j;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwLin_AdjLevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwLin_AdjLevelSet(void);
  
  /*!
   * \brief Compute the upwind flux between two nodes i and j.
   * \param[out] val_residual_i - Pointer to the total residual at node i.
   * \param[out] val_residual_j - Pointer to the total residual at node j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method from node i to node j (implicit computation).
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method from node j to node i (implicit computation).
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                       su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CUpwLin_AdjTurb
 * \brief Class for performing a linear upwind solver for the adjoint turbulence equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CUpwLin_AdjTurb : public CNumerics {
private:
  su2double *Velocity_i;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwLin_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwLin_AdjTurb(void);
  
  /*!
   * \brief Compute the adjoint upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_TurbSA
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaral turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CUpwSca_TurbSA : public CNumerics {
private:
  su2double *Velocity_i, *Velocity_j;
  bool implicit, grid_movement, incompressible;
  su2double q_ij, a0, a1;
  unsigned short iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TurbSA(void);
  
  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_TurbML
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaral turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CUpwSca_TurbML : public CNumerics {
private:
  su2double *Velocity_i, *Velocity_j;
  bool implicit, grid_movement, incompressible;
  su2double q_ij, a0, a1;
  unsigned short iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbML(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TurbML(void);
  
  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_TurbSST
 * \brief Class for doing a scalar upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 * \version 4.0.1 "Cardinal"
 */
class CUpwSca_TurbSST : public CNumerics {
private:
  su2double *Velocity_i, *Velocity_j;
  bool implicit, grid_movement, incompressible;
  su2double Density_i, Density_j,
  q_ij,
  a0, a1;
  unsigned short iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TurbSST(void);
  
  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_TransLM
 * \brief Class for doing a scalar upwind solver for the Spalart-Allmaras turbulence model equations with transition.
 * \ingroup ConvDiscr
 * \author A. Aranake.
 * \version 4.0.1 "Cardinal"
 */
class CUpwSca_TransLM : public CNumerics {
private:
  su2double *Velocity_i, *Velocity_j;
  bool implicit, grid_movement;
  su2double q_ij, a0, a1;
  unsigned short iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TransLM(void);
  
  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_AdjTurb
 * \brief Class for doing a scalar upwind solver for the adjoint turbulence equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CUpwSca_AdjTurb : public CNumerics {
private:
  su2double *Velocity_i, *Velocity_j;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_AdjTurb(void);
  
  /*!
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};


/*!
 * \class CCentJST_Flow
 * \brief Class for centered shceme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentJST_KE_Flow : public CNumerics {
  
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Epsilon_2, cte_0, cte_1, /*!< \brief Artificial dissipation values. */
  ProjGridVel;  /*!< \brief Projected grid velocity. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_KE_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_KE_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentJST_Flow
 * \brief Class for centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentJST_Flow : public CNumerics {
  
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Epsilon_2, Epsilon_4, cte_0, cte_1, /*!< \brief Artificial dissipation values. */
  ProjGridVel;  /*!< \brief Projected grid velocity. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentJSTArtComp_Flow
 * \brief Class for centered scheme - JST (artificial compressibility).
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentJSTArtComp_Flow : public CNumerics {
  
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanBetaInc2, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Epsilon_2, Epsilon_4, cte_0, cte_1; /*!< \brief Artificial dissipation values. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement, /*!< \brief Modification for grid movement. */
  gravity; /*!< \brief computation with gravity force. */
  su2double Froude; /*!< \brief Froude number. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJSTArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJSTArtComp_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentJST_AdjFlow
 * \brief Class for and adjoint centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentJST_AdjFlow : public CNumerics {
private:
  su2double *Diff_Psi, *Diff_Lapl;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanPhi;
  unsigned short iDim, jDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j, ProjPhi, ProjPhi_Vel, sq_vel, phis1, phis2;
  su2double MeanPsiRho, MeanPsiE, Param_p, Param_Kappa_4, Param_Kappa_2, Local_Lambda_i, Local_Lambda_j, MeanLambda;
  su2double Phi_i, Phi_j, sc4, StretchingFactor, Epsilon_4, Epsilon_2;
  bool implicit, grid_movement;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a JST method.
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CCentJSTArtComp_AdjFlow
 * \brief Class for and adjoint centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentJSTArtComp_AdjFlow : public CNumerics {
private:
  su2double sc2, *Diff_Psi, *Diff_Lapl;
  su2double *Velocity_i, *Velocity_j;
  su2double **Proj_Jac_Tensor_i, **Proj_Jac_Tensor_j;
  unsigned short iDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j;
  su2double Param_p, Param_Kappa_4, Param_Kappa_2, Local_Lambda_i, Local_Lambda_j, MeanLambda;
  su2double Phi_i, Phi_j, sc4, StretchingFactor, Epsilon_4, Epsilon_2;
  bool implicit;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJSTArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJSTArtComp_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a JST method.
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CCentJST_LinFlow
 * \brief Class for linearized centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentJST_LinFlow : public CNumerics {
private:
  su2double *Diff_DeltaU, *Diff_Lapl;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanDeltaVel, *MeanVelocity;
  su2double **MeanJacobian;
  su2double **Jacobian_i, **Jacobian_j;
  unsigned short iDim, iVar, jVar;
  su2double sq_vel, Density_i, DensityEnergy_i, Energy_i, Pressure_i, Density_j, DensityEnergy_j, Energy_j,
  Pressure_j, Param_p, Param_Kappa_4, Local_Lambda_i, Local_Lambda_j, MeanLambda, sc4, StretchingFactor,
  Epsilon_4, MeanDeltaRho, MeanDeltaE, ProjVelocity_i, ProjVelocity_j, MeanDensity, MeanPressure,
  MeanEnthalpy, MeanEnergy, Phi_i, Phi_j;
  
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_LinFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_LinFlow(void);
  
  /*!
   * \brief Compute the linearized flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_resconv, su2double *val_resvisc, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                        CConfig *config);
};

/*!
 * \class CCentLax_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentLax_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, /*!< \brief Difference of conservative variables. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_0, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc0, StretchingFactor, /*!< \brief Streching parameters. */
  Epsilon_0, cte; /*!< \brief Artificial dissipation values. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  su2double ProjGridVel;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLax_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentLaxArtComp_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme (artificial compressibility).
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentLaxArtComp_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, /*!< \brief Difference of conservative variables. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanBetaInc2, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_0, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc0, StretchingFactor, /*!< \brief Streching parameters. */
  Epsilon_0; /*!< \brief Artificial dissipation values. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement, /*!< \brief Modification for grid movement. */
  gravity; /*!< \brief Modification for for gravity force. */
  su2double Froude;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLaxArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLaxArtComp_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentLax_AdjFlow
 * \brief Class for computing the Lax-Friedrich adjoint centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentLax_AdjFlow : public CNumerics {
private:
  su2double *Diff_Psi;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanPhi;
  unsigned short iDim, jDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j, ProjPhi, ProjPhi_Vel, sq_vel, phis1, phis2,
  MeanPsiRho, MeanPsiE, Param_p, Param_Kappa_0, Local_Lambda_i, Local_Lambda_j, MeanLambda,
  Phi_i, Phi_j, sc2, StretchingFactor, Epsilon_0;
  bool implicit, grid_movement;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLax_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLax_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a Lax method.
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CCentLaxArtComp_AdjFlow
 * \brief Class for computing the Lax-Friedrich adjoint centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentLaxArtComp_AdjFlow : public CNumerics {
private:
  su2double *Diff_Psi;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanPhi, **Proj_Jac_Tensor_i, **Proj_Jac_Tensor_j;
  unsigned short iDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j, Param_p, Param_Kappa_0,
  Local_Lambda_i, Local_Lambda_j, MeanLambda,
  Phi_i, Phi_j, sc2, StretchingFactor, Epsilon_0;
  bool implicit;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLaxArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLaxArtComp_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a Lax method.
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CCentLax_LinFlow
 * \brief Class for computing the Lax-Friedrich linearized centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CCentLax_LinFlow : public CNumerics {
private:
  su2double *Diff_DeltaU;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanDeltaVel, *MeanVelocity;
  su2double **MeanJacobian;
  su2double **Jacobian_i;
  su2double **Jacobian_j;
  unsigned short iDim, iVar, jVar;
  su2double sq_vel, Density_i, DensityEnergy_i, Energy_i, Pressure_i, Density_j,
  DensityEnergy_j, Energy_j, Pressure_j, Param_p, Param_Kappa_0,
  Local_Lambda_i, Local_Lambda_j, MeanLambda, StretchingFactor,
  Epsilon_i, MeanDeltaRho, MeanDeltaE, ProjVelocity_i, ProjVelocity_j,
  MeanDensity, MeanPressure,
  MeanEnthalpy, MeanEnergy, Phi_i, Phi_j,
  sc2;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLax_LinFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLax_LinFlow(void);
  
  /*!
   * \brief Compute the linearized flow residual using a Lax method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_resconv, su2double *val_resvisc, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;	   /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,				   /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,				   /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,					   /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,                /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,                   /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,				/*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij;						/*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

class CAvgGrad_Flow2 : public CNumerics {
private:
	unsigned short iDim, iVar, jVar;	   /*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar,				   /*!< \brief Mean primitive variables. */
	*PrimVar_i, *PrimVar_j,				   /*!< \brief Primitives variables at point i and 1. */
	Total_Viscosity_i,Total_Viscosity_j,
	  **Mean_GradPrimVar,**Mean_Tau,visc_div_vel,					   /*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity,                /*!< \brief Mean value of the viscosity. */
	Mean_Eddy_Viscosity,                   /*!< \brief Mean value of the eddy viscosity. */
	Mean_Thermal_Conductivity,             /*!< \brief Mean value of the thermal conductivity. */
	Mean_Cp,                               /*!< \brief Mean value of the Cp. */
	Mean_turb_ke,				/*!< \brief Mean value of the turbulent kinetic energy. */
	*ProjFlux,	/*!< \brief Projection of the viscous fluxes. */
	dist_ij;						/*!< \brief Length of the edge and face. */
	bool implicit; /*!< \brief Implicit calculus. */

public:
    
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_Flow2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGrad_Flow2(void);
    
	/*!
	 * \brief Compute the viscous flow residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CGeneralAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author M.Pini, S. Vitale
 * \version 3.2.1 "eagle"
 */

class CGeneralAvgGrad_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;	   /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,				   /*!< \brief Mean primitive variables. */
  *Mean_SecVar,        				   /*!< \brief Mean secondary variables. */
  *PrimVar_i, *PrimVar_j,				   /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,					   /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,                /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,                   /*!< \brief Mean value of the eddy viscosity. */
  Mean_Thermal_Conductivity,             /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp,                               /*!< \brief Mean value of the Cp. */
  Mean_turb_ke,				/*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij;						/*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralAvgGrad_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradArtComp_Flow
 * \brief Class for computing viscous term using an average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradArtComp_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;	/*!< \brief Iterators in dimension an variable. */
  su2double **Mean_GradPrimVar,					/*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
  dist_ij;							/*!< \brief Length of the edge and face. */
  bool implicit;				/*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradArtComp_Flow(void);
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbSA
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_TurbSA : public CNumerics {
private:
  
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;
  su2double nu_i, nu_j, nu_e;
  su2double dist_ij_2;
  su2double proj_vector_ij;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSA(void);
  
  /*!
   * \brief Compute the viscous turbulence terms residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbSA_Neg
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_TurbSA_Neg : public CNumerics {
private:
  
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;
  su2double cn1, fn, Xi;
  su2double nu_i, nu_j, nu_ij, nu_tilde_ij, nu_e;
  su2double dist_ij_2;
  su2double proj_vector_ij;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSA_Neg(void);
  
  /*!
   * \brief Compute the viscous turbulence terms residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbML
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_TurbML : public CNumerics {
private:
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;
  su2double nu_i, nu_j, nu_e;
  su2double dist_ij_2;
  su2double proj_vector_ij;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbML(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbML(void);
  
  /*!
   * \brief Compute the viscous turbulence terms residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TransLM
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_TransLM : public CNumerics {
private:
  su2double **Mean_GradTransVar;
  su2double *Proj_Mean_GradTransVar_Kappa, *Proj_Mean_GradTransVar_Edge;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;
  su2double dist_ij_2;
  su2double proj_vector_ij;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TransLM(void);
  
  /*!
   * \brief Compute the viscous turbulence terms residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_AdjFlow
 * \brief Class for computing the adjoint viscous terms.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_AdjFlow : public CNumerics {
private:
  su2double *Velocity_i;	/*!< \brief Auxiliary vector for storing the velocity of point i. */
  su2double *Velocity_j;	/*!< \brief Auxiliary vector for storing the velocity of point j. */
  su2double *Mean_Velocity;
  su2double *Mean_GradPsiE;	/*!< \brief Counter for dimensions of the problem. */
  su2double **Mean_GradPhi;	/*!< \brief Counter for dimensions of the problem. */
  su2double *Edge_Vector;	/*!< \brief Vector going from node i to node j. */
  bool implicit;			/*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_AdjFlow(void);
  
  /*!
   * \brief Residual computation.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGradArtComp_AdjFlow
 * \brief Class for computing the adjoint viscous terms.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradArtComp_AdjFlow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;	/*!< \brief Iterators in dimension an variable. */
  su2double **Mean_GradPsiVar,					/*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
  dist_ij;							/*!< \brief Length of the edge and face. */
  bool implicit;				/*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradArtComp_AdjFlow(void);
  
  /*!
   * \brief Residual computation.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;		/*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,				/*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,									/*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,	/*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,      /*!< \brief Mean value of the laminar viscosity. */
  Mean_Eddy_Viscosity,         /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,				 /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij_2;					 /*!< \brief Length of the edge and face. */
  bool implicit;			/*!< \brief Implicit calculus. */
  bool limiter;			/*!< \brief Viscous limiter. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 * \version 3.2.9 "eagle"
 */
class CAvgGradCorrected_Flow2 : public CNumerics {
private:
	unsigned short iDim, iVar, jVar;		/*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
	*PrimVar_i, *PrimVar_j,				/*!< \brief Primitives variables at point i and 1. */
	Total_Viscosity_i, Total_Viscosity_j,   
	*Edge_Vector,									/*!< \brief Vector form point i to point j. */
	  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,**Mean_Tau,visc_div_vel,	/*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity,      /*!< \brief Mean value of the laminar viscosity. */
	Mean_Eddy_Viscosity,         /*!< \brief Mean value of the eddy viscosity. */
	Mean_Thermal_Conductivity,   /*!< \brief Mean value of the thermal conductivity. */
	Mean_Cp,                     /*!< \brief Mean value of the specific heat. */
	Mean_turb_ke,				 /*!< \brief Mean value of the turbulent kinetic energy. */
	dist_ij_2,					 /*!< \brief Length of the edge and face. */
	*ProjFlux;	/*!< \brief Projection of the viscous fluxes. */
	bool implicit;			/*!< \brief Implicit calculus. */
  bool limiter;			/*!< \brief Viscous limiter. */

public:
    
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradCorrected_Flow2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
    
	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGradCorrected_Flow2(void);
    
	/*!
	 * \brief Compute the viscous flow residual using an average of gradients with correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};


/*!
 * \class CGeneralAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author M. Pini, S. Vitale
 * \version 3.2.1 "eagle"
 */
class CGeneralAvgGradCorrected_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;		/*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
  *Mean_SecVar,			        		/*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,			    	/*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,									/*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,	/*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,      /*!< \brief Mean value of the laminar viscosity. */
  Mean_Eddy_Viscosity,         /*!< \brief Mean value of the eddy viscosity. */
  Mean_Thermal_Conductivity,   /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp,                     /*!< \brief Mean value of the specific heat. */
  Mean_turb_ke,				 /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij_2;					 /*!< \brief Length of the edge and face. */
  bool implicit;			/*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralAvgGradCorrected_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrectedArtComp_Flow
 * \brief Class for computing viscous term using an average of gradients with correction (artificial compresibility).
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrectedArtComp_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;	/*!< \brief Iterators in dimension an variable. */
  su2double *PrimVar_i, *PrimVar_j,			/*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,								/*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,	/*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,			/*!< \brief Mean value of the viscosity. */
  dist_ij_2;					/*!< \brief Length of the edge and face. */
  bool implicit;			/*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrectedArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrectedArtComp_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_TurbSA
 * \brief Class for computing viscous term using average of gradients with correction (Spalart-Allmaras turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_TurbSA : public CNumerics {
private:
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma, nu_i, nu_j, nu_e, dist_ij_2, proj_vector_ij;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_TurbSA(void);
  
  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_TurbSA_Neg
 * \brief Class for computing viscous term using average of gradients with correction (Spalart-Allmaras turbulence model).
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_TurbSA_Neg : public CNumerics {
private:
  
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
  su2double *Edge_Vector;
  su2double sigma;
  su2double cn1, fn, Xi;
  su2double nu_ij, nu_tilde_ij;
  bool implicit, incompressible;
  su2double nu_i, nu_j, nu_e, dist_ij_2, proj_vector_ij;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_TurbSA_Neg(void);
  
  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_TurbML
 * \brief Class for computing viscous term using average of gradients with correction (Spalart-Allmaras turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_TurbML : public CNumerics {
private:
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma, nu_i, nu_j, nu_e, dist_ij_2, proj_vector_ij;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_TurbML(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_TurbML(void);
  
  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_TransLM
 * \brief Class for computing viscous term using average of gradients with correction (Spalart-Allmaras turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_TransLM : public CNumerics {
private:
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_TransLM(void);
  
  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_TurbSST : public CNumerics {
private:
  su2double sigma_k1,                     /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
  sigma_k2,
  sigma_om1,
  sigma_om2;
  
  su2double diff_kine,                     /*!< \brief Diffusivity for viscous terms of tke eq */
  diff_omega;                           /*!< \brief Diffusivity for viscous terms of omega eq */
  
  su2double *Edge_Vector,                  /*!< \brief Vector from node i to node j. */
  dist_ij_2,                            /*!< \brief |Edge_Vector|^2 */
  proj_vector_ij;                       /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  
  su2double **Mean_GradTurbVar,            /*!< \brief Average of gradients at cell face */
  *Proj_Mean_GradTurbVar_Normal,        /*!< \brief Mean_gradTurbVar DOT normal */
  *Proj_Mean_GradTurbVar_Edge,          /*!< \brief Mean_gradTurbVar DOT Edge_Vector */
  *Proj_Mean_GradTurbVar_Corrected;
  
  su2double F1_i, F1_j;                    /*!< \brief Menter's first blending function */
  
  bool implicit, incompressible;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double* constants, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSST(void);
  
  /*!
   * \brief Sets value of first blending function.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j) { F1_i = val_F1_i; F1_j = val_F1_j;}
  
  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients wtih correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
  
};

/*!
 * \class CAvgGradCorrected_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_TurbSST : public CNumerics {
private:
  su2double sigma_k1,                     /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
  sigma_k2,
  sigma_om1,
  sigma_om2;
  
  su2double diff_kine,                     /*!< \brief Diffusivity for viscous terms of tke eq */
  diff_omega;                           /*!< \brief Diffusivity for viscous terms of omega eq */
  
  su2double *Edge_Vector,                  /*!< \brief Vector from node i to node j. */
  dist_ij_2,                            /*!< \brief |Edge_Vector|^2 */
  proj_vector_ij;                       /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  
  su2double **Mean_GradTurbVar,            /*!< \brief Average of gradients at cell face */
  *Proj_Mean_GradTurbVar_Normal,        /*!< \brief Mean_gradTurbVar DOT normal */
  *Proj_Mean_GradTurbVar_Edge,          /*!< \brief Mean_gradTurbVar DOT Edge_Vector */
  *Proj_Mean_GradTurbVar_Corrected;
  
  su2double F1_i, F1_j;                    /*!< \brief Menter's first blending function */
  
  bool implicit, incompressible;
  unsigned short iVar, iDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double* constants, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_TurbSST(void);
  
  /*!
   * \brief Sets value of first blending function.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j) { F1_i = val_F1_i; F1_j = val_F1_j;}
  
  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients wtih correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
  
};

/*!
 * \class CAvgGradCorrected_AdjFlow
 * \brief Class for computing the adjoint viscous terms, including correction.
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_AdjFlow : public CNumerics {
private:
  su2double *Velocity_i;	/*!< \brief Auxiliary vector for storing the velocity of point i. */
  su2double *Velocity_j;	/*!< \brief Auxiliary vector for storing the velocity of point j. */
  su2double *Mean_Velocity;
  su2double **Mean_GradPsiVar;	/*!< \brief Counter for dimensions of the problem. */
  su2double *Edge_Vector;	/*!< \brief Vector going from node i to node j. */
  su2double *Proj_Mean_GradPsiVar_Edge;	/*!< \brief Projection of Mean_GradPsiVar onto Edge_Vector. */
  su2double *Mean_GradPsiE;	/*!< \brief Counter for dimensions of the problem. */
  su2double **Mean_GradPhi;	/*!< \brief Counter for dimensions of the problem. */
  bool implicit;          /*!< \brief Boolean controlling Jacobian calculations. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow viscous residual in a non-conservative way using an average of gradients and derivative correction.
   * \param[out] val_residual_i - Pointer to the viscous residual at point i.
   * \param[out] val_residual_j - Pointer to the viscous residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGradCorrectedArtComp_AdjFlow
 * \brief Class for computing the adjoint viscous terms, including correction.
 * \ingroup ViscDiscr
 * \author F.Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrectedArtComp_AdjFlow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;	/*!< \brief Iterators in dimension an variable. */
  su2double *PsiVar_i, *PsiVar_j,			/*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,								/*!< \brief Vector form point i to point j. */
  **Mean_GradPsiVar, *Proj_Mean_GradPsiVar_Edge,	/*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,			/*!< \brief Mean value of the viscosity. */
  dist_ij_2;					/*!< \brief Length of the edge and face. */
  bool implicit;			/*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrectedArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrectedArtComp_AdjFlow(void);
  
  /*!
   * \brief Compute the adjoint flow viscous residual in a non-conservative way using an average of gradients and derivative correction.
   * \param[out] val_residual_i - Pointer to the viscous residual at point i.
   * \param[out] val_residual_j - Pointer to the viscous residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                        su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_AdjTurb
 * \brief Class for adjoint turbulent using average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_AdjTurb : public CNumerics {
private:
  su2double **Mean_GradTurbPsi;
  su2double *Proj_Mean_GradTurbPsi_Kappa, *Proj_Mean_GradTurbPsi_Edge, *Proj_Mean_GradTurbPsi_Corrected;
  su2double *Edge_Vector;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_AdjTurb(void);
  
  /*!
   * \brief Compute the adjoint turbulent residual using average of gradients and a derivative correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CAvgGrad_AdjTurb
 * \brief Class for adjoint turbulent using average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_AdjTurb : public CNumerics {
private:
  su2double **Mean_GradTurbPsi;
  su2double *Proj_Mean_GradTurbPsi_Kappa, *Proj_Mean_GradTurbPsi_Edge, *Proj_Mean_GradTurbPsi_Corrected;
  su2double *Edge_Vector;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_AdjTurb(void);
  
  /*!
   * \brief Compute the adjoint turbulent residual using average of gradients and a derivative correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CGalerkin_Flow
 * \brief Class for computing the stiffness matrix of the Galerkin method.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CGalerkin_Flow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGalerkin_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CGalerkin_Flow(void);
  
  /*!
   * \brief Computing stiffness matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double **val_stiffmatrix_elem, CConfig *config);
};

/*!
 * \class CGalerkin_FEA
 * \brief Class for computing the stiffness matrix of the Galerkin method.
 * \ingroup ViscDiscr
 * \author F. Palacios, R.Sanchez
 * \version 4.0.1 "Cardinal"
 */
class CGalerkin_FEA : public CNumerics {
  
  su2double E;				/*!< \brief Young's modulus of elasticity. */
  su2double Nu;			/*!< \brief Poisson's ratio. */
  su2double Rho_s;		/*!< \brief Structural density. */
  su2double Mu;			/*!< \brief Lame's coeficient. */
  su2double Lambda;	/*!< \brief Lame's coeficient. */

public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGalerkin_FEA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CGalerkin_FEA(void);
  
  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Fnodal - Forces at the nodes in cartesian coordinates.
   * \param[in] Pnodal - Pressure at the nodes.
   * \param[in] CoordCorners[2][2] - Coordiantes of the corners.
   */
  void PressInt_Linear(su2double CoordCorners[4][3], su2double *tn_e, su2double Fnodal[12]);
  
  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Tau_0 - Stress tensor at the node 0.
   * \param[in] Tau_1 - Stress tensor at the node 1.
   * \param[in] Fnodal - Forces at the nodes in cartesian coordinates.
   * \param[in] CoordCorners[2][2] - Coordiantes of the corners.
   */
  void ViscTermInt_Linear(su2double CoordCorners[2][2], su2double Tau_0[3][3], su2double Tau_1[3][3],  su2double FviscNodal[4]);
  
  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
   * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
   */
  su2double ShapeFunc_Triangle(su2double Xi, su2double Eta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
   * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
   */
  su2double ShapeFunc_Quadrilateral(su2double Xi, su2double Eta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
   * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
   */
  su2double ShapeFunc_Tetra(su2double Xi, su2double Eta, su2double Mu, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
   * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
   */
  su2double ShapeFunc_Prism(su2double Xi, su2double Eta, su2double Mu, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
   * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
   */
  su2double ShapeFunc_Pyram(su2double Xi, su2double Eta, su2double Mu, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
   * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
   * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
   */
  su2double ShapeFunc_Hexa(su2double Xi, su2double Eta, su2double Mu, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]);
  
  /*!
   * \brief Computing stiffness matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFEA_StiffMatrix2D(su2double **StiffMatrix_Elem, su2double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d);
  
  /*!
   * \brief Computing stiffness matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFEA_StiffMatrix3D(su2double **StiffMatrix_Elem, su2double CoordCorners[8][3], unsigned short nNodes);
  
  /*!
   * \brief Computing mass matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFEA_StiffMassMatrix2D(su2double **StiffMatrix_Elem, su2double **MassMatrix_Elem, su2double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d);
  
  /*!
   * \brief Computing mass matrix of the Galerkin method.
   * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFEA_StiffMassMatrix3D(su2double **StiffMatrix_Elem, su2double **MassMatrix_Elem, su2double CoordCorners[8][3], unsigned short nNodes);
  
  /*!
   * \brief Computing stresses in FEA method at the nodes.
   * \param[in] config - Definition of the particular problem.
   */
  void GetFEA_StressNodal2D(su2double StressVector[8][3], su2double DispElement[8], su2double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d);
  
  
  /*!
   * \brief Computing stresses in FEA method at the nodes.
   * \param[in] config - Definition of the particular problem.
   */
  void GetFEA_StressNodal3D(su2double StressVector[8][6], su2double DispElement[24], su2double CoordCorners[8][3], unsigned short nNodes);
  
  /*!
   * \brief Computing dead load vector of the Galerkin method.
   * \param[out] val_deadloadvector_elem - Dead load at the nodes for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFEA_DeadLoad2D(su2double *DeadLoadVector_Elem, su2double CoordCorners[8][3], unsigned short nNodes, su2double matDensity);
  
  /*!
   * \brief Computing stiffness matrix of the Galerkin method.
   * \param[out] val_deadloadvector_elem - Dead load at the nodes for Galerkin computation.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFEA_DeadLoad3D(su2double *DeadLoadVector_Elem, su2double CoordCorners[8][3], unsigned short nNodes, su2double matDensity);
  
};

/*!
 * \class CSourceNothing
 * \brief Dummy class.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CSourceNothing : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceNothing(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceNothing(void);
};

/*!
 * \class CSourcePieceWise_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CSourcePieceWise_TurbSA : public CNumerics {
private:
  su2double cv1_3;
  su2double k2;
  su2double cb1;
  su2double cw2;
  su2double ct3;
  su2double ct4;
  su2double cw3_6;
  su2double cb2_sigma;
  su2double sigma;
  su2double cb2;
  su2double cw1;
  unsigned short iDim;
  su2double nu, Ji, fv1, fv2, ft2, Omega, S, Shat, inv_Shat, dist_i_2, Ji_2, Ji_3, inv_k2_d2;
  su2double r, g, g_6, glim, fw;
  su2double norm2_Grad;
  su2double dfv1, dfv2, dShat;
  su2double dr, dg, dfw;
  bool incompressible;
  bool rotating_frame;
  su2double intermittency;
  su2double Production, Destruction, CrossProduction;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TurbSA(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] intermittency_in - Value of the intermittency.
   */
  void SetIntermittency(su2double intermittency_in);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  void SetProduction(su2double val_production);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  void SetDestruction(su2double val_destruction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  void SetCrossProduction(su2double val_crossproduction);
  
  /*!
   * \brief ______________.
   */
  su2double GetProduction(void);
  
  /*!
   * \brief  ______________.
   */
  su2double GetDestruction(void);
  
  /*!
   * \brief  ______________.
   */
  su2double GetCrossProduction(void);
};

/*!
 * \class CSourcePieceWise_TurbSA_Neg
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CSourcePieceWise_TurbSA_Neg : public CNumerics {
private:
  su2double cv1_3;
  su2double k2;
  su2double cb1;
  su2double cw2;
  su2double ct3;
  su2double ct4;
  su2double cw3_6;
  su2double cb2_sigma;
  su2double sigma;
  su2double cb2;
  su2double cw1;
  unsigned short iDim;
  su2double nu, Ji, fv1, fv2, ft2, Omega, S, Shat, inv_Shat, dist_i_2, Ji_2, Ji_3, inv_k2_d2;
  su2double r, g, g_6, glim, fw;
  su2double norm2_Grad;
  su2double dfv1, dfv2, dShat;
  su2double dr, dg, dfw;
  bool incompressible;
  bool rotating_frame;
  su2double intermittency;
  su2double Production, Destruction, CrossProduction;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TurbSA_Neg(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] intermittency_in - Value of the intermittency.
   */
  void SetIntermittency(su2double intermittency_in);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  void SetProduction(su2double val_production);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  void SetDestruction(su2double val_destruction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  void SetCrossProduction(su2double val_crossproduction);
  
  /*!
   * \brief ______________.
   */
  su2double GetProduction(void);
  
  /*!
   * \brief  ______________.
   */
  su2double GetDestruction(void);
  
  /*!
   * \brief  ______________.
   */
  su2double GetCrossProduction(void);
};

/*!
 * \class CSourcePieceWise_TurbML
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CSourcePieceWise_TurbML : public CNumerics {
private:
  bool incompressible;
  bool transition;
  bool rotating_frame;
  su2double intermittency;
  su2double Production, Destruction, CrossProduction;
  CScalePredictor* MLModel;
  
  su2double uInfinity;
  
  SpalartAllmarasInputs* SAInputs;
  SpalartAllmarasConstants* SAConstants;
  
  int nResidual;
  int nJacobian;
  
  string featureset;

  su2double** DUiDXj;
  su2double* DNuhatDXj;
  
public:
  
  bool isInBL;
  su2double fw;
  su2double fWake;
  SpalartAllmarasOtherOutputs* SAOtherOutputs;
  
  su2double *SAResidual;
  su2double * SANondimResidual;
  su2double* Residual;
  su2double * NondimResidual;
  su2double *ResidualDiff;
  su2double *NondimResidualDiff;
  su2double* SAJacobian;
  CSANondimInputs* SANondimInputs;
  su2double NuhatGradNorm;
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbML(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TurbML(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] intermittency_in - Value of the intermittency.
   */
  void SetIntermittency(su2double intermittency_in);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  void SetProduction(su2double val_production);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  void SetDestruction(su2double val_destruction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  void SetCrossProduction(su2double val_crossproduction);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  su2double GetProduction(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  su2double GetDestruction(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  su2double GetCrossProduction(void);
  
  su2double SAProduction, SADestruction, SACrossProduction, SASource, MLProduction, MLDestruction, MLCrossProduction, MLSource, SourceDiff;
  
  int NumResidual();
};

/*!
 * \class CSourcePieceWise_TransLM
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CSourcePieceWise_TransLM : public CNumerics {
private:
  
  /*-- SA model constants --*/
  su2double cv1_3;
  su2double k2;
  su2double cb1;
  su2double cw2;
  su2double cw3_6;
  su2double sigma;
  su2double cb2;
  su2double cw1;
  
  /*-- gamma-theta model constants --*/
  su2double c_e1;
  su2double c_a1;
  su2double c_e2;
  su2double c_a2;
  su2double sigmaf;
  su2double s1;
  su2double c_theta;
  su2double sigmat;
  
  /*-- Correlation constants --*/
  su2double flen_global;
  su2double alpha_global;
  su2double Vorticity;

  bool implicit;
  
public:
  bool debugme; // For debugging only, remove this. -AA
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TransLM(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual_TransLM(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config, su2double &gamma_sep);
  
  void CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(su2double *TransVar_i, su2double *TransVar_id, su2double *val_residual, su2double *val_residuald, CConfig *config);
};

/*!
 * \class CSourcePieceWise_TurbSST
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author A. Campos.
 * \version 4.0.1 "Cardinal"
 */
class CSourcePieceWise_TurbSST : public CNumerics {
private:
  su2double F1_i,
  F1_j,
  F2_i,
  F2_j;
  
  su2double alfa_1,
  alfa_2,
  beta_1,
  beta_2,
  sigma_omega_1,
  sigma_omega_2,
  beta_star,
  a1;
  
  su2double CDkw_i, CDkw_j;
  
  bool incompressible;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double* constants, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TurbSST(void);
  
  /*!
   * \brief Set the value of the first blending function.
   * \param[in] val_F1_i - Value of the first blending function at point i.
   * \param[in] val_F1_j - Value of the first blending function at point j.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j);
  
  /*!
   * \brief Set the value of the second blending function.
   * \param[in] val_F2_i - Value of the second blending function at point i.
   * \param[in] val_F2_j - Value of the second blending function at point j.
   */
  void SetF2blending(su2double val_F2_i, su2double val_F2_j);
  
  /*!
   * \brief Set the value of the cross diffusion for the SST model.
   * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
   * \param[in] val_CDkw_j - Value of the cross diffusion at point j.
   */
  virtual void SetCrossDiff(su2double val_CDkw_i, su2double val_CDkw_j);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
};

/*!
 * \class CSourceGravity
 * \brief Class for the source term integration of the gravity force.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CSourceGravity : public CNumerics {
  su2double Froude;
  bool compressible, incompressible, freesurface;
  
public:
  
  /*!
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceGravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceGravity(void);
  
  /*!
   * \brief Source term integration for the poissonal potential.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config);
};

/*!
 * \class CSourceViscous_AdjFlow
 * \brief Class for source term integration in adjoint problem.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CSourceViscous_AdjFlow : public CNumerics {
private:
  su2double *Velocity, *GradDensity, *GradInvDensity, *dPoDensity2, *alpha, *beta, *Sigma_5_vec;
  su2double **GradVel_o_Rho, **sigma, **Sigma_phi, **Sigma_5_Tensor, **Sigma;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceViscous_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceViscous_AdjFlow(void);
  
  /*!
   * \brief Source term integration of the flow adjoint equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, CConfig *config);
  
};

/*!
 * \class CSourcePieceWise_AdjTurb
 * \brief Class for source term integration of the adjoint turbulent equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CSourcePieceWise_AdjTurb : public CNumerics {
private:
  su2double **tau, *Velocity;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_AdjTurb(void);
  
  /*!
   * \brief Source term integration of the adjoint turbulence equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CSourcePieceWise_AdjLevelSet
 * \brief Class for source term integration of the adjoint level set equation.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CSourcePieceWise_AdjLevelSet : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_AdjLevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_AdjLevelSet(void);
  
  /*!
   * \brief Source term integration of the adjoint poisson potential equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config);
};

/*!
 * \class CSourceConservative_AdjFlow
 * \brief Class for source term integration in adjoint problem using a conservative scheme.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CSourceConservative_AdjFlow : public CNumerics {
private:
  su2double *Velocity, *Residual_i, *Residual_j, *Mean_Residual;
  su2double **Mean_PrimVar_Grad;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceConservative_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceConservative_AdjFlow(void);
  
  /*!
   * \brief Source term integration using a conservative scheme.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config);
};

/*!
 * \class CSourceConservative_AdjTurb
 * \brief Class for source term integration in adjoint turbulent problem using a conservative scheme.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 4.0.1 "Cardinal"
 */
class CSourceConservative_AdjTurb : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceConservative_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceConservative_AdjTurb(void);
  
  /*!
   * \brief Source term integration using a conservative scheme.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CSourceRotatingFrame_Flow
 * \brief Class for a rotating frame source term.
 * \ingroup SourceDiscr
 * \author F. Palacios, T. Economon.
 * \version 4.0.1 "Cardinal"
 */
class CSourceRotatingFrame_Flow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceRotatingFrame_Flow(void);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
};

/*!
 * \class CSourceRotatingFrame_AdjFlow
 * \brief Source term class for rotating frame adjoint.
 * \ingroup SourceDiscr
 * \author T. Economon.
 * \version 4.0.1 "Cardinal"
 */
class CSourceRotatingFrame_AdjFlow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceRotatingFrame_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceRotatingFrame_AdjFlow(void);
  
  /*!
   * \brief Residual of the adjoint rotating frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
};

/*!
 * \class CSourceAxisymmetric_Flow
 * \brief Class for source term for solving axisymmetric problems.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CSourceAxisymmetric_Flow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceAxisymmetric_Flow(void);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config);
  
};

/*!
 * \class CSourceAxisymmetric_AdjFlow
 * \brief Class for source term for solving axisymmetric problems.
 * \ingroup SourceDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CSourceAxisymmetric_AdjFlow : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceAxisymmetric_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceAxisymmetric_AdjFlow(void);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config);
};

/*!
 * \class CSourceWindGust
 * \brief Class for a source term due to a wind gust.
 * \ingroup SourceDiscr
 * \author S. Padrón
 * \version 4.0.1 "Cardinal"
 */
class CSourceWindGust : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceWindGust(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceWindGust(void);
  
  /*!
   * \brief Residual of the wind gust source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
};

/*!
 * \class CSource_Template
 * \brief Dummy class.
 * \ingroup SourceDiscr
 * \author A. Lonkar.
 * \version 4.0.1 "Cardinal"
 */
class CSource_Template : public CNumerics {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config -  Name of the input config file
   *
   */
  CSource_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSource_Template(void);
};

/*!
 * \class CConvectiveTemplate
 * \brief Class for setting up new method for spatial discretization of convective terms in flow Equations
 * \ingroup ConvDiscr
 * \author A. Lonkar
 * \version 4.0.1 "Cardinal"
 */
class CConvective_Template : public CNumerics {
private:
  
  /* define private variables here */
  bool implicit;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
  unsigned short iDim, iVar, jVar, kVar;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CConvective_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CConvective_Template(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CViscous_Template
 * \brief Class for computing viscous term using average of gradients.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CViscous_Template : public CNumerics {
private:
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CViscous_Template(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CViscous_Template(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwRoe_TNE2
 * \brief Class for evaluating the Riemann problem using Roe's scheme for a two-temperature model.
 * \ingroup ConvDiscr
 * \author S. R. Copeland
 * \version 2.0.6
 */
class CUpwRoe_TNE2 : public CNumerics {
private:
  bool implicit, ionization;
  su2double *Diff_U;
  su2double *RoeU, *RoeV;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double RoeSoundSpeed;
  su2double ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  su2double Proj_ModJac_Tensor_ij, R;
  su2double *RoedPdU;
  unsigned short nSpecies, nPrimVar, nPrimVarGrad, nVar, nDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoe_TNE2(unsigned short val_nDim, unsigned short val_nVar,
               unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
               CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_TNE2(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
};


/*!
 * \class CUpwMSW_TNE2
 * \brief Class for solving a flux-vector splitting method by Steger & Warming, modified version.
 * \ingroup ConvDiscr
 * \author S. Copeland
 * \version 4.0.1 "Cardinal"
 */
class CUpwMSW_TNE2 : public CNumerics {
private:
  bool ionization, implicit;
  su2double *Diff_U;
  su2double *u_i, *u_j, *ust_i, *ust_j;
  su2double *Fc_i, *Fc_j;
  su2double *Lambda_i, *Lambda_j;
  su2double *rhos_i, *rhos_j, *rhosst_i, *rhosst_j;
  su2double *Ust_i, *Ust_j, *Vst_i, *Vst_j;
  su2double *dPdUst_i, *dPdUst_j;
  su2double **P_Tensor, **invP_Tensor;
  unsigned short nSpecies, nPrimVar, nPrimVarGrad, nVar, nDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar
   * \param[in] val_nPrimVarGrad
   * \param[in] config - Definition of the particular problem.
   */
  CUpwMSW_TNE2(unsigned short val_nDim, unsigned short val_nVar,
               unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
               CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwMSW_TNE2(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
};

/*!
 * \class CUpwAUSM_TNE2
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 2.0.6
 */
class CUpwAUSM_TNE2 : public CNumerics {
private:
  bool implicit, ionization;
  su2double *FcL, *FcR, *FcLR;
  su2double *dmLP, *dmRM, *dpLP, *dpRM;
  su2double *daL, *daR;
  su2double *rhos_i, *u_i;
  su2double *rhos_j, *u_j;
  su2double a_i, P_i, h_i, ProjVel_i;
  su2double a_j, P_j, h_j, ProjVel_j;
  unsigned short nSpecies, nVar, nDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSM_TNE2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSM_TNE2(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwAUSM_TNE2
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 2.0.6
 */
class CUpwAUSMPWplus_TNE2 : public CNumerics {
private:
  bool implicit, ionization;
  su2double *FcL, *FcR;
  su2double *dmLdL, *dmLdR, *dmRdL, *dmRdR;
  su2double *dmLPdL, *dmLPdR, *dmRMdL, *dmRMdR;
  su2double *dmbLPdL, *dmbLPdR, *dmbRMdL, *dmbRMdR;
  su2double *dpLPdL, *dpLPdR, *dpRMdL, *dpRMdR;
  su2double *dHnL, *dHnR;
  su2double *daL, *daR;
  su2double *rhos_i, *u_i;
  su2double *rhos_j, *u_j;
  su2double *dPdU_i, *dPdU_j;
  unsigned short nSpecies, nVar, nDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPWplus_TNE2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSMPWplus_TNE2(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};


/*!
 * \class CCentLax_TNE2
 * \brief Class for computing the Lax-Friedrich centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 2.0.6
 */
class CCentLax_TNE2 : public CNumerics {
private:
  unsigned short jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U; /*!< \brief Difference of conservative variables. */
  su2double *MeanU, *MeanV;
  su2double *MeandPdU;
  su2double *ProjFlux;  /*!< \brief Projected inviscid flux tensor. */
  su2double Param_p, Param_Kappa_0; /*!< \brief Artificial dissipation parameters. */
  su2double Local_Lambda_i, Local_Lambda_j, MeanLambda; /*!< \brief Local eigenvalues. */
  su2double Phi_i, Phi_j, sc0, StretchingFactor; /*!< \brief Streching parameters. */
  su2double Epsilon_0, cte; /*!< \brief Artificial dissipation values. */
  bool implicit; /*!< \brief Implicit time integration. */
  bool ionization;  /*!< \brief Charged species with the mixture. */
  unsigned short nSpecies, nVar, nPrimVar, nPrimVarGrad, nDim;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLax_TNE2(unsigned short val_nDim, unsigned short val_nVar,
                unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLax_TNE2(void);
  
  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_resconv, su2double *val_resvisc, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author S. R. Copeland
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_TNE2 : public CNumerics {
private:
  unsigned short iDim, iVar, nPrimVar, nPrimVarGrad;		/*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,				/*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,						/*!< \brief Mean value of the gradient. */
  *Mean_Diffusion_Coeff, /*!< \brief Mean value of the species diffusion coefficient. */
  Mean_Laminar_Viscosity, /*!< \brief Mean value of the viscosity. */
  Mean_Thermal_Conductivity, /*!< \brief Mean value of the thermal conductivity. */
  Mean_Thermal_Conductivity_ve, /*!< \brief Mean value of the vib-el. thermal conductivity. */
  dist_ij;						/*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of variables in the primitive variable gradient.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TNE2(unsigned short val_nDim,
                unsigned short val_nVar,
                unsigned short val_nPrimVar,
                unsigned short val_nPrimVarGrad,
                CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TNE2(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual,
                       su2double **val_Jacobian_i,
                       su2double **val_Jacobian_j,
                       CConfig *config);
};


/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author S. R. Copeland
 * \version 4.0.1 "Cardinal"
 */
class CAvgGradCorrected_TNE2 : public CNumerics {
private:
  unsigned short iDim, iVar, nPrimVar, nPrimVarGrad;		/*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,				/*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,						/*!< \brief Mean value of the gradient. */
  *Edge_Vector,
  *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  *Mean_Diffusion_Coeff, /*!< \brief Mean value of the species diffusion coefficient. */
  Mean_Laminar_Viscosity, /*!< \brief Mean value of the viscosity. */
  Mean_Thermal_Conductivity, /*!< \brief Mean value of the thermal conductivity. */
  Mean_Thermal_Conductivity_ve, /*!< \brief Mean value of the vib-el. thermal conductivity. */
    dist_ij;						/*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of variables in the primitive variable gradient.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_TNE2(unsigned short val_nDim,
                         unsigned short val_nVar,
                         unsigned short val_nPrimVar,
                         unsigned short val_nPrimVarGrad,
                         CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_TNE2(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual,
                       su2double **val_Jacobian_i,
                       su2double **val_Jacobian_j,
                       CConfig *config);
};


/*!
 * \class CSource_TNE2
 * \brief Class for two-temperature model source terms.
 * \ingroup SourceDiscr
 * \author S. Copeland
 * \version 2.0.6
 */
class CSource_TNE2 : public CNumerics {
private:
  bool   implicit, ionization;
  unsigned short nSpecies, nVar, nPrimVar, nPrimVarGrad;
  int    *alphak, *betak;
  su2double *X; // Mole fraction
  su2double **RxnConstantTable;
  su2double *dkf, *dkb, *dRfok, *dRbok, *A;
  su2double *eves, *Cvvs, *Cves;
  //  CVariable *var;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSource_TNE2(unsigned short val_nDim,
               unsigned short val_nVar,
               unsigned short val_nPrimVar,
               unsigned short val_nPrimVarGrad,
               CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSource_TNE2(void);
  
  /*!
   * \brief Source residual of the chemistry.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeChemistry(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
  
  /*!
   * \brief Calculates constants used for Keq correlation.
   * \param[out] A - Pointer to coefficient array.
   * \param[in] val_reaction - Reaction number indicator.
   * \param[in] config - Definition of the particular problem.
   */
  void GetKeqConstants(su2double *A, unsigned short val_reaction, CConfig *config);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeVibRelaxation(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config);
};


/*!
 * \class CUpwRoe_AdjTNE2
 * \brief Class for solving an approximate Riemann solver of Roe
 *        for the adjoint flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 2.0.6
 */
class CUpwRoe_AdjTNE2 : public CNumerics {
private:
  bool implicit;
  unsigned short nVar, nPrimVar, nPrimVarGrad, nSpecies;
  su2double *MeanU, *MeanV, *MeandPdU;
  su2double *DiffPsi;
  su2double *UnitNormal;
  su2double *Lambda;
  su2double **Ai, **Aj;
  su2double **P, **invP, **PLPinv;
  
  //  CVariable *var;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoe_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar,
                  unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                  CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_AdjTNE2(void);
  
  /*!
   * \brief Compute the adjoint Roe's flux between two nodes i and j.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                       CConfig *config);
};

/*!
 * \class CUpwSW_AdjTNE2
 * \brief Class for solving an approximate Riemann solver of Roe
 *        for the adjoint flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 2.0.6
 */
class CUpwSW_AdjTNE2 : public CNumerics {
private:
  bool implicit;
  unsigned short nVar, nPrimVar, nPrimVarGrad, nSpecies;
  su2double *DiffPsi;
  su2double *UnitNormal;
  su2double *Lambda_i, *Lambda_j;
  su2double **P, **invP, **PLPinv;
  su2double **Ai, **Aj;
  
  //  CVariable *var;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSW_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar,
                 unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                 CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSW_AdjTNE2(void);
  
  /*!
   * \brief Compute the adjoint Roe's flux between two nodes i and j.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                       CConfig *config);
};


/*!
 * \class CCentJST_AdjTNE2
 * \brief Class for and adjoint centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 2.0.6
 */
class CCentJST_AdjTNE2 : public CNumerics {
private:
  su2double *Diff_Psi, *Diff_Lapl;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanPhi;
  unsigned short iDim, jDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j, ProjPhi, ProjPhi_Vel, sq_vel, phis1, phis2;
  su2double MeanPsiRho, MeanPsiE, Param_p, Param_Kappa_4, Param_Kappa_2, Local_Lambda_i, Local_Lambda_j, MeanLambda;
  su2double Phi_i, Phi_j, sc4, StretchingFactor, Epsilon_4, Epsilon_2;
  bool implicit, grid_movement, rotating_frame;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_AdjTNE2(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a JST method.
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};


/*!
 * \class CCentLax_AdjTNE2
 * \brief Class for computing the Lax-Friedrich adjoint centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 2.0.6
 */
class CCentLax_AdjTNE2 : public CNumerics {
private:
  bool implicit;
  unsigned short nVar, nPrimVar, nPrimVarGrad, nSpecies, nDim;
  su2double *DiffPsi, *MeanPsi;
  su2double Param_p, Param_Kappa_0;
  su2double **Proj_Jac_Tensor_i, **Proj_Jac_Tensor_j;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLax_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar,
                   unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                   CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLax_AdjTNE2(void);
  
  /*!
   * \brief Compute the adjoint flow residual using a Lax method.
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i,
                        su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                        su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};

/*!
 * \class CAvgGrad_AdjTNE2
 * \brief Class for computing the adjoint viscous terms.
 * \ingroup ViscDiscr
 * \author F. Palacios
 * \version 4.0.1 "Cardinal"
 */
class CAvgGrad_AdjTNE2 : public CNumerics {
private:
  su2double *vel, *vel_i, *vel_j;
  su2double *Mean_GradPsiE;	/*!< \brief Mean gradient in the adjoint  energy between nodes i and j. */
  su2double *Mean_GradPsiEve; /*!< \brief Mean gradient in the adjoint vibrational energy between nodes i and j. */
  su2double **Mean_GradPhi;	/*!< \brief Counter for dimensions of the problem. */
  su2double *Edge_Vector;	/*!< \brief Vector going from node i to node j. */
  su2double **SigmaPhi;
  su2double **SigmaPsiE;
  bool implicit;			/*!< \brief Implicit calculus. */
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_AdjTNE2(void);
  
  /*!
   * \brief Residual computation.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};

/*!
 * \class CSource_AdjTNE2
 * \brief Class for adjoint two-temperature model source terms.
 * \ingroup SourceDiscr
 * \author S. Copeland
 * \version 2.0.6
 */
class CSource_AdjTNE2 : public CNumerics {
private:
  bool   implicit;
  unsigned short nSpecies, nVar, nPrimVar, nPrimVarGrad;
  su2double *rhos, *vel;
  su2double *GInvRho, **GVeloRho, **tau, **eta, **pi, **zeta;
  su2double *GPhiGInvRho, *GPsiEZetaTau;
  su2double **Av2, **Av3, **Av4;
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSource_AdjTNE2(unsigned short val_nDim,
                  unsigned short val_nVar,
                  unsigned short val_nPrimVar,
                  unsigned short val_nPrimVarGrad,
                  CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSource_AdjTNE2(void);
  
  /*!
   * \brief Source residual of the chemistry.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeSourceViscous(su2double *val_residual, CConfig *config);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeSourceConservative(su2double *val_residual, CConfig *config);
};


#include "numerics_structure.inl"
