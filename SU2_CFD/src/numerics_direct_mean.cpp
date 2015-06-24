/*!
 * \file numerics_direct_mean.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios, T. Economon
 * \version 4.0.0 "Cardinal"
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

#include "../include/numerics_structure.hpp"
#include <limits>

CCentJST_Flow::CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  grid_movement = config->GetGrid_Movement();
  
  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new double [nVar];
  Diff_Lapl = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  ProjFlux = new double [nVar];
  
}

CCentJST_Flow::~CCentJST_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
}

void CCentJST_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                    CConfig *config) {
  
  double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;
  
  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }

  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);
  
  /*--- Residual of the inviscid flux ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  
  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Adjustment due to grid motion ---*/
  
  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar] + U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }
  
  /*--- Computes differences btw. Laplacians and conservative variables,
   with a correction for the enthalpy ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  }
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  /*--- Compute the local spectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  /*--- Adjustment due to mesh motion ---*/
  
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
  sc4 = sc2*sc2/4.0;
  
  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/
  
  if (implicit) {
    
    cte_0 = (Epsilon_2 + Epsilon_4*double(Neighbor_i+1))*StretchingFactor*MeanLambda;
    cte_1 = (Epsilon_2 + Epsilon_4*double(Neighbor_j+1))*StretchingFactor*MeanLambda;
    
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_1;
    }
    
    /*--- Last row of Jacobian_i ---*/
    
    val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;
    
    /*--- Last row of Jacobian_j ---*/
    
    val_Jacobian_j[nVar-1][0] -= cte_1*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte_1*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte_1*Gamma;
    
  }
  
}

CCentJST_KE_Flow::CCentJST_KE_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  grid_movement = config->GetGrid_Movement();

  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();

  /*--- Allocate some structures ---*/
  Diff_U = new double [nVar];
  Diff_Lapl = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  ProjFlux = new double [nVar];

}

CCentJST_KE_Flow::~CCentJST_KE_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
}

void CCentJST_KE_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                    CConfig *config) {

  double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/

  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;

  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }

  /*--- Recompute conservative variables ---*/

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;

  /*--- Compute mean values of the variables ---*/

  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);

  /*--- Get projected flux tensor ---*/

  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);

  /*--- Residual of the inviscid flux ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];

  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }

  /*--- Adjustment due to grid motion ---*/

  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }

  /*--- Computes differences btw. Laplacians and conservative variables,
   with a correction for the enthalpy ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  }
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;

  /*--- Compute the local spectral radius and the stretching factor ---*/

  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*--- Adjustment due to mesh motion ---*/

  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);

  sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
  sc4 = sc2*sc2/4.0;

  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;

  /*--- Compute viscous part of the residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
      val_residual[iVar] += Epsilon_2*(Diff_U[iVar])*StretchingFactor*MeanLambda;

  /*--- Jacobian computation ---*/

  if (implicit) {

    cte_0 = Epsilon_2*StretchingFactor*MeanLambda;

    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_0;
    }

    /*--- Last row of Jacobian_i ---*/

    val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;

    /*--- Last row of Jacobian_j ---*/

    val_Jacobian_j[nVar-1][0] -= cte_1*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte_1*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte_1*Gamma;

  }

}


CCentLax_Flow::CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  
  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  ProjFlux = new double [nVar];
  
}

CCentLax_Flow::~CCentLax_Flow(void) {
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
  
}

void CCentLax_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                    CConfig *config) {
  
  double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;
  
  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  
  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);
  
  /*--- Residual of the inviscid flux ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  
  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Adjustment due to grid motion ---*/
  
  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }
  
  /*--- Computes differences btw. conservative variables,
   with a correction for the enthalpy ---*/
  
  for (iVar = 0; iVar < nDim+1; iVar++)
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  Diff_U[nDim+1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  /*--- Compute the local spectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  /*--- Adjustment due to grid motion ---*/
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc0 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*double(nDim)/3.0;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/
  
  if (implicit) {
    cte = Epsilon_0*StretchingFactor*MeanLambda;
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte;
      val_Jacobian_j[iVar][iVar] -= cte;
    }
    
    /*--- Last row of Jacobian_i ---*/
    
    val_Jacobian_i[nVar-1][0] += cte*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte*Gamma;
    
    /*--- Last row of Jacobian_j ---*/
    
    val_Jacobian_j[nVar-1][0] -= cte*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte*Gamma;
    
  }
  
}

CUpwCUSP_Flow::CUpwCUSP_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  grid_movement = config->GetGrid_Movement();
  
  /*--- Allocate some structures ---*/
  Diff_U = new double [nVar];
  Diff_Flux = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  ProjFlux = new double [nVar];
  ProjFlux_i = new double [nVar];
  ProjFlux_j = new double [nVar];
  Jacobian = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian[iVar] = new double [nVar];
  }
}

CUpwCUSP_Flow::~CUpwCUSP_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Flux;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian[iVar];
  }
}

void CUpwCUSP_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                     CConfig *config) {
  
  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;

  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  
  SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel_i));
  SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel_j));

  /*-- Face area ---*/
  
  Area = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	/*-- Unit normal ---*/
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
    ProjVelocity +=  MeanVelocity[iDim]*UnitNormal[iDim];
  }
  MeanSoundSpeed = 0.5*(SoundSpeed_i+SoundSpeed_j);
  MeanEnergy = 0.5*(Energy_i+Energy_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);
  
  /*--- Residual of the inviscid flux ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  
  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Computes differences conservative variables,
   with a correction for the enthalpy ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  /*--- Computes differences projected fluxes,
   with a correction for the enthalpy ---*/
  
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, UnitNormal, ProjFlux_i);
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, UnitNormal, ProjFlux_j);
  
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_Flux[iVar] = ProjFlux_i[iVar]-ProjFlux_j[iVar];
  
  /*--- Compute dissipation parameters ---*/
  
  Mach = ProjVelocity / MeanSoundSpeed;
  
  LamdaNeg = ProjVelocity - MeanSoundSpeed;
  LamdaPos = ProjVelocity + MeanSoundSpeed;
  
  if ((0.0 <= Mach) && (Mach < 1.0)) Beta = + max(0.0, (ProjVelocity + LamdaNeg)/(ProjVelocity - LamdaNeg));
  if ((-1.0 <= Mach) && (Mach < 0.0)) Beta = - max(0.0, (ProjVelocity + LamdaPos)/(ProjVelocity - LamdaPos));
  if (fabs(Mach) >= 1.0) Beta = Mach/fabs(Mach);
  
  if (Beta == 0.0) Nu_c = fabs(ProjVelocity);
  if ((Beta > 0.0) && ((0.0 < Mach) && (Mach < 1.0))) Nu_c = - (1.0-Beta)*LamdaNeg;
  if ((Beta < 0.0) && ((-1.0 < Mach) && (Mach < 0.0))) Nu_c = (1.0-Beta)*LamdaPos;
  if (fabs(Mach) >= 1) Nu_c = 0.0;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += (0.5*Nu_c*Diff_U[iVar] + 0.5*Beta*Diff_Flux[iVar])*Area;

  /*--- Jacobian computation ---*/
  
  if (implicit) {
    
    cte_0 = 0.5*Nu_c*Area;
    cte_1 = 0.5*Beta*Area;
    
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_0;
    }
    
    /*--- Last row of Jacobian_i (solution difference contribution) ---*/
    
    val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;
    
    /*--- Last row of Jacobian_j (solution difference contribution) ---*/
    
    val_Jacobian_j[nVar-1][0] -= cte_0*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte_0*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte_0*Gamma;
    
    /*--- Flux difference contribution ---*/
    
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 1.0, Jacobian);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] += cte_1*Jacobian[iVar][jVar];
    
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 1.0, Jacobian);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] -= cte_1*Jacobian[iVar][jVar];
    
  }
  
}

CUpwAUSM_Flow::CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  delta_vel  = new double [nDim];
  delta_wave = new double [nVar];
  ProjFlux_i = new double [nVar];
  ProjFlux_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
}

CUpwAUSM_Flow::~CUpwAUSM_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwAUSM_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel)));
  
  /*--- Primitive variables at point j ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel)));
  
  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  mL	= ProjVelocity_i/SoundSpeed_i;
  mR	= ProjVelocity_j/SoundSpeed_j;
  
  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0);
  else mLP = 0.5*(mL+fabs(mL));
  
  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0);
  else mRM = 0.5*(mR-fabs(mR));
  
  mF = mLP + mRM;
  
  if (fabs(mL) <= 1.0) pLP = 0.25*Pressure_i*(mL+1.0)*(mL+1.0)*(2.0-mL);
  else pLP = 0.5*Pressure_i*(mL+fabs(mL))/mL;
  
  if (fabs(mR) <= 1.0) pRM = 0.25*Pressure_j*(mR-1.0)*(mR-1.0)*(2.0+mR);
  else pRM = 0.5*Pressure_j*(mR-fabs(mR))/mR;
  
  pF = pLP + pRM;
  Phi = fabs(mF);
  
  val_residual[0] = 0.5*(mF*((Density_i*SoundSpeed_i)+(Density_j*SoundSpeed_j))-Phi*((Density_j*SoundSpeed_j)-(Density_i*SoundSpeed_i)));
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = 0.5*(mF*((Density_i*SoundSpeed_i*Velocity_i[iDim])+(Density_j*SoundSpeed_j*Velocity_j[iDim]))
                                -Phi*((Density_j*SoundSpeed_j*Velocity_j[iDim])-(Density_i*SoundSpeed_i*Velocity_i[iDim])))+UnitNormal[iDim]*pF;
  val_residual[nVar-1] = 0.5*(mF*((Density_i*SoundSpeed_i*Enthalpy_i)+(Density_j*SoundSpeed_j*Enthalpy_j))-Phi*((Density_j*SoundSpeed_j*Enthalpy_j)-(Density_i*SoundSpeed_i*Enthalpy_i)));
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  /*--- Roe's Jacobian for AUSM (this must be fixed) ---*/
  if (implicit) {
    
    /*--- Mean Roe variables iPoint and jPoint ---*/
    R = sqrt(fabs(Density_j/Density_i));
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
}

CUpwHLLC_Flow::CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  delta_vel  = new double [nDim];
  delta_wave = new double [nVar];
  ProjFlux_i = new double [nVar];
  ProjFlux_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
  
}

CUpwHLLC_Flow::~CUpwHLLC_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwHLLC_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  sq_vel_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel_i += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel_i)));
  
  /*--- Primitive variables at point j ---*/
  
  sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_j += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel_j)));
  
  /*--- Projected velocities ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Roe's aveaging ---*/
  
  Rrho = sqrt(fabs(Density_j/Density_i));
  tmp = 1.0/(1.0+Rrho);
  for (iDim = 0; iDim < nDim; iDim++)
    velRoe[iDim] = tmp*(Velocity_i[iDim] + Velocity_j[iDim]*Rrho);
  
  uRoe  = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    uRoe += velRoe[iDim]*UnitNormal[iDim];
  
  gamPdivRho = tmp*((Gamma*Pressure_i/Density_i+0.5*(Gamma-1.0)*sq_vel_i) + (Gamma*Pressure_j/Density_j+0.5*(Gamma-1.0)*sq_vel_j)*Rrho);
  sq_velRoe = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sq_velRoe += velRoe[iDim]*velRoe[iDim];
  
  cRoe  = sqrt(fabs(gamPdivRho - ((Gamma+Gamma)*0.5-1.0)*0.5*sq_velRoe));
  
  /*--- Speed of sound at L and R ---*/
  
  sL = min(uRoe-cRoe, ProjVelocity_i-SoundSpeed_i);
  sR = max(uRoe+cRoe, ProjVelocity_j+SoundSpeed_j);
  
  /*--- speed of contact surface ---*/
  
  sM = (Pressure_i-Pressure_j
        - Density_i*ProjVelocity_i*(sL-ProjVelocity_i)
        + Density_j*ProjVelocity_j*(sR-ProjVelocity_j))
  /(Density_j*(sR-ProjVelocity_j)-Density_i*(sL-ProjVelocity_i));
  
  /*--- Pressure at right and left (Pressure_j=Pressure_i) side of contact surface ---*/
  
  pStar = Density_j * (ProjVelocity_j-sR)*(ProjVelocity_j-sM) + Pressure_j;
  
  if (sM >= 0.0) {
    if (sL > 0.0) {
      val_residual[0] = Density_i*ProjVelocity_i;
      for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = Density_i*Velocity_i[iDim]*ProjVelocity_i + Pressure_i*UnitNormal[iDim];
      val_residual[nVar-1] = Energy_i*Density_i*ProjVelocity_i + Pressure_i*ProjVelocity_i;
    }
    else {
      invSLmSs = 1.0/(sL-sM);
      sLmuL = sL-ProjVelocity_i;
      rhoSL = Density_i*sLmuL*invSLmSs;
      for (iDim = 0; iDim < nDim; iDim++)
        rhouSL[iDim] = (Density_i*Velocity_i[iDim]*sLmuL+(pStar-Pressure_i)*UnitNormal[iDim])*invSLmSs;
      eSL = (sLmuL*Energy_i*Density_i-Pressure_i*ProjVelocity_i+pStar*sM)*invSLmSs;
      
      val_residual[0] = rhoSL*sM;
      for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = rhouSL[iDim]*sM + pStar*UnitNormal[iDim];
      val_residual[nVar-1] = (eSL+pStar)*sM;
    }
  }
  else {
    if (sR >= 0.0) {
      invSRmSs = 1.0/(sR-sM);
      sRmuR = sR-ProjVelocity_j;
      rhoSR = Density_j*sRmuR*invSRmSs;
      for (iDim = 0; iDim < nDim; iDim++)
        rhouSR[iDim] = (Density_j*Velocity_j[iDim]*sRmuR+(pStar-Pressure_j)*UnitNormal[iDim])*invSRmSs;
      eSR = (sRmuR*Energy_j*Density_j-Pressure_j*ProjVelocity_j+pStar*sM)*invSRmSs;
      
      val_residual[0] = rhoSR*sM;
      for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = rhouSR[iDim]*sM + pStar*UnitNormal[iDim];
      val_residual[nVar-1] = (eSR+pStar)*sM;
    }
    else {
      val_residual[0] = Density_j*ProjVelocity_j;
      for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = Density_j*Velocity_j[iDim]*ProjVelocity_j + Pressure_j*UnitNormal[iDim];
      val_residual[nVar-1] = Energy_j*Density_j*ProjVelocity_j + Pressure_j*ProjVelocity_j;
    }
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  if (implicit) {
    
    /*--- Mean Roe variables iPoint and jPoint ---*/
    
    R = sqrt(fabs(Density_j/Density_i));
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
  
}

#ifdef CHECK

int UgpWithCvCompFlow::calcEulerFluxMatrices_HLLC(double (*val_Jacobian_i)[5], double (*val_Jacobian_j)[5], double (*val_Jacobian_i_Scal)[6], double (*val_Jacobian_j_Scal)[6],
                                                  const double Density_i, const double *uL, const double pL, const double TL, const double h0, const double RL, const double gammaL, const double *scalL, const double kL,
                                                  const double Density_j, const double *uR, const double pR, const double TR, const double h1, const double RR, const double gammaR, const double *scalR, const double kR,
                                                  const double area, const double *nVec, const int nScal, const double surfVeloc)
{

  double unL  = vecDotVec3d(uL, nVec);
  double uLuL = vecDotVec3d(uL, uL);
  double cL   = sqrt(gammaL*pL/Density_i);
  double hL   = gammaL/(gammaL-1.0)*pL/Density_i + 0.5*uLuL + kL;
  //  double hL   = h0 + 0.5*uLuL + kL;
  double eL   = hL*Density_i-pL;
  
  double unR  = vecDotVec3d(uR, nVec);
  double uRuR = vecDotVec3d(uR, uR);
  double cR   = sqrt(gammaR*pR/Density_j);
  double hR   = gammaR/(gammaR-1.0)*pR/Density_j + 0.5*uRuR + kR;
  //  double hR   = h1 + 0.5*uRuR + kR;
  double eR   = hR*Density_j-pR;
  
  
  // Roe's aveaging
  double Rrho = sqrt(Density_j/Density_i);
  double tmp = 1.0/(1.0+Rrho);
  double velRoe[3];
  for (int i=0; i<3; i++)
    velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
  double uRoe  = vecDotVec3d(velRoe, nVec);
  double hRoe = tmp*(hL + hR*Rrho);
  
  //  double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d(velRoe, velRoe)));
  double gamPdivRho = tmp*((gammaL*pL/Density_i+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/Density_j+0.5*(gammaR-1.0)*uRuR)*Rrho);
  double cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*vecDotVec3d(velRoe, velRoe));
  
  // speed of sound at L and R
  double sL = min(uRoe-cRoe, unL-cL);
  double sR = max(uRoe+cRoe, unR+cR);
  
  // speed of contact surface
  double sM = (pL-pR-Density_i*unL*(sL-unL)+Density_j*unR*(sR-unR))/(Density_j*(sR-unR)-Density_i*(sL-unL));
  
  // pressure at right and left (pR=pL) side of contact surface
  double pStar = Density_j*(unR-sR)*(unR-sM)+pR;
  
  if (sM >= 0.0) {
    
    if (sL > 0.0) {
      
      double nVecArea[3];
      for (int i=0; i<3; i++) nVecArea[i] = nVec[i]*area;
      
      calcJacobianA(val_Jacobian_i, uL, pL, Density_i, nVecArea, 0.5*(gammaL+gammaR), 0.0);
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] = 0.0;
      
    }
    else {
      
      double invSLmSs = 1.0/(sL-sM);
      double sLmuL = sL-unL;
      double rhoSL = Density_i*sLmuL*invSLmSs;
      double rhouSL[3];
      
      for (int i=0; i<3; i++)
        rhouSL[i] = (Density_i*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
      
      double eSL = (sLmuL*eL-pL*unL+pStar*sM)*invSLmSs;
      double gammaLM1 = (gammaL-1.0);
      double gammaRM1 = (gammaR-1.0);
      double invrhotld = 1.0/(Density_j*(sR-unR)-Density_i*(sL-unL));
      
      double dSMdUL[5], dSMdUR[5];
      double dpsdUL[5], dpsdUR[5];
      
      dSMdUL[0] = -unL*unL + uLuL*gammaLM1/2.0 + sM*sL;
      dSMdUL[1] =  nVec[0]*(2.0*unL-sL-sM) - gammaLM1*uL[0];
      dSMdUL[2] =  nVec[1]*(2.0*unL-sL-sM) - gammaLM1*uL[1];
      dSMdUL[3] =  nVec[2]*(2.0*unL-sL-sM) - gammaLM1*uL[2];
      dSMdUL[4] =  gammaLM1;
      
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dSMdUL[i] *= invrhotld;
        dpsdUL[i] = Density_j*(sR-unR)*dSMdUL[i];
      }
      
      dSMdUR[0] =  unR*unR - uRuR*gammaRM1/2.0 - sM*sR;
      dSMdUR[1] = -nVec[0]*(2.0*unR-sR-sM) + gammaRM1*uR[0];
      dSMdUR[2] = -nVec[1]*(2.0*unR-sR-sM) + gammaRM1*uR[1];
      dSMdUR[3] = -nVec[2]*(2.0*unR-sR-sM) + gammaRM1*uR[2];
      dSMdUR[4] = -gammaRM1;
      
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dSMdUR[i] *= invrhotld;
        dpsdUR[i] = Density_i*(sL-unL)*dSMdUR[i];
      }
      
      calcSubSonicJacobeanHLLC(val_Jacobian_i, val_Jacobian_j,
                               Density_i, uL, pL, eL, unL, uLuL, sL,
                               rhoSL, rhouSL, eSL, dSMdUL,
                               dSMdUR, dpsdUL, dpsdUR, sM, pStar, 0.5*(gammaL+gammaR), nVec);
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[0][i] =  val_Jacobian_i[0][i]*sM + dSMdUL[i]*rhoSL;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[1][i] =  val_Jacobian_i[1][i]*sM + dSMdUL[i]*rhouSL[0] + dpsdUL[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[2][i] =  val_Jacobian_i[2][i]*sM + dSMdUL[i]*rhouSL[1] + dpsdUL[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[3][i] =  val_Jacobian_i[3][i]*sM + dSMdUL[i]*rhouSL[2] + dpsdUL[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[4][i] = (val_Jacobian_i[4][i]+dpsdUL[i])*sM + (eSL+pStar)*dSMdUL[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] *= area;
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[0][i] =  val_Jacobian_j[0][i]*sM + dSMdUR[i]*rhoSL;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[1][i] =  val_Jacobian_j[1][i]*sM + dSMdUR[i]*rhouSL[0] + dpsdUR[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[2][i] =  val_Jacobian_j[2][i]*sM + dSMdUR[i]*rhouSL[1] + dpsdUR[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[3][i] =  val_Jacobian_j[3][i]*sM + dSMdUR[i]*rhouSL[2] + dpsdUR[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[4][i] = (val_Jacobian_j[4][i]+dpsdUR[i])*sM + (eSL+pStar)*dSMdUR[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] *= area;
      
    }
  }
  
  else {
    
    if (sR >= 0.0) {
      
      double invSRmSs = 1.0/(sR-sM);
      double sRmuR = sR-unR;
      double rhoSR = Density_j*sRmuR*invSRmSs;
      double rhouSR[3];
      for (int i=0; i<3; i++)
        rhouSR[i] = (Density_j*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
      double eSR = (sRmuR*eR-pR*unR+pStar*sM)*invSRmSs;
      double gammaLM1 = (gammaL-1.0);
      double gammaRM1 = (gammaR-1.0);
      double invrhotld = 1.0/(Density_j*(sR-unR)-Density_i*(sL-unL));
      
      double dSMdUL[5], dSMdUR[5];
      double dpsdUL[5], dpsdUR[5];
      
      dSMdUL[0] = -unL*unL + uLuL*gammaLM1/2.0 + sM*sL;
      dSMdUL[1] =  nVec[0]*(2.0*unL-sL-sM) - gammaLM1*uL[0];
      dSMdUL[2] =  nVec[1]*(2.0*unL-sL-sM) - gammaLM1*uL[1];
      dSMdUL[3] =  nVec[2]*(2.0*unL-sL-sM) - gammaLM1*uL[2];
      dSMdUL[4] =  gammaLM1;
      
      for (iVar = 0; iVar < nVar; iVar++) {
        dSMdUL[i] *= invrhotld;
        dpsdUL[i] = Density_j*(sR-unR)*dSMdUL[i];
      }
      
      dSMdUR[0] =  unR*unR - uRuR*gammaRM1/2.0 - sM*sR;
      dSMdUR[1] = -nVec[0]*(2.0*unR-sR-sM) + gammaRM1*uR[0];
      dSMdUR[2] = -nVec[1]*(2.0*unR-sR-sM) + gammaRM1*uR[1];
      dSMdUR[3] = -nVec[2]*(2.0*unR-sR-sM) + gammaRM1*uR[2];
      dSMdUR[4] = -gammaRM1;
      
      for (iVar = 0; iVar < nVar; iVar++) {
        dSMdUR[i] *= invrhotld;
        dpsdUR[i] = Density_i*(sL-unL)*dSMdUR[i];
      }
      
      calcSubSonicJacobeanHLLC(val_Jacobian_j, val_Jacobian_i,
                               Density_j, uR, pR, eR, unR, uRuR, sR,
                               rhoSR, rhouSR, eSR,
                               dSMdUR, dSMdUL, dpsdUR, dpsdUL, sM, pStar, 0.5*(gammaL+gammaR), nVec);
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[0][i] =  val_Jacobian_i[0][i]*sM + dSMdUL[i]*rhoSR;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[1][i] =  val_Jacobian_i[1][i]*sM + dSMdUL[i]*rhouSR[0] + dpsdUL[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[2][i] =  val_Jacobian_i[2][i]*sM + dSMdUL[i]*rhouSR[1] + dpsdUL[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[3][i] =  val_Jacobian_i[3][i]*sM + dSMdUL[i]*rhouSR[2] + dpsdUL[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[4][i] = (val_Jacobian_i[4][i]+dpsdUL[i])*sM + (eSR+pStar)*dSMdUL[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] *= area;
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[0][i] =  val_Jacobian_j[0][i]*sM + dSMdUR[i]*rhoSR;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[1][i] =  val_Jacobian_j[1][i]*sM + dSMdUR[i]*rhouSR[0] + dpsdUR[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[2][i] =  val_Jacobian_j[2][i]*sM + dSMdUR[i]*rhouSR[1] + dpsdUR[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[3][i] =  val_Jacobian_j[3][i]*sM + dSMdUR[i]*rhouSR[2] + dpsdUR[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[4][i] = (val_Jacobian_j[4][i]+dpsdUR[i])*sM + (eSR+pStar)*dSMdUR[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] *= area;
      
    }
    
    else {
      
      double nVecArea[3];
      for (int i=0; i<3; i++)        nVecArea[i] = nVec[i]*area;
      calcJacobianA(val_Jacobian_j, uR, pR, Density_j, nVecArea, 0.5*(gammaL+gammaR), 0.0);
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] = 0.0;
      
    }
    
  }
  
}

void UgpWithCvCompFlow::calcSubSonicJacobeanHLLC(double (*AL)[5], double (*AR)[5],
                                                 double Density_i, const double *uL, double pL, double eL, double qL, double psiL, double SL,
                                                 double rhoSL, double *rhouSL, double eSL,
                                                 double *dSMdUL, double *dSMdUR, double *dpsdUL, double *dpsdUR, double SM, double pS,
                                                 double gamma, const double *nV) // nV is not normalized
{
  
  double gammaMinus1 = (gamma-1.0);
  double omL = 1.0/(SL-SM);
  
  AL[0][0] =  SL    + rhoSL*dSMdUL[0];
  AL[0][1] = -nV[0] + rhoSL*dSMdUL[1];
  AL[0][2] = -nV[1] + rhoSL*dSMdUL[2];
  AL[0][3] = -nV[2] + rhoSL*dSMdUL[3];
  AL[0][4] =        + rhoSL*dSMdUL[4];
  
  AL[1][0] =    qL*uL[0]       - nV[0]*psiL*gammaMinus1/2.0   + nV[0]*dpsdUL[0] + rhouSL[0]*dSMdUL[0];
  AL[1][1] =  SL - qL          + nV[0]*(gamma-2.0)*uL[0]      + nV[0]*dpsdUL[1] + rhouSL[0]*dSMdUL[1];
  AL[1][2] =     - uL[0]*nV[1] + nV[0]*gammaMinus1*uL[1]      + nV[0]*dpsdUL[2] + rhouSL[0]*dSMdUL[2];
  AL[1][3] =     - uL[0]*nV[2] + nV[0]*gammaMinus1*uL[2]      + nV[0]*dpsdUL[3] + rhouSL[0]*dSMdUL[3];
  AL[1][4] = -gammaMinus1*nV[0]                               + nV[0]*dpsdUL[4] + rhouSL[0]*dSMdUL[4];
  
  AL[2][0] =    qL*uL[1]       - nV[1]*psiL*gammaMinus1/2.0   + nV[1]*dpsdUL[0] + rhouSL[1]*dSMdUL[0];
  AL[2][1] =     - uL[1]*nV[0] + nV[1]*gammaMinus1*uL[0]      + nV[1]*dpsdUL[1] + rhouSL[1]*dSMdUL[1];
  AL[2][2] =  SL - qL          + nV[1]*(gamma-2.0)*uL[1]      + nV[1]*dpsdUL[2] + rhouSL[1]*dSMdUL[2];
  AL[2][3] =     - uL[1]*nV[2] + nV[1]*gammaMinus1*uL[2]      + nV[1]*dpsdUL[3] + rhouSL[1]*dSMdUL[3];
  AL[2][4] = -gammaMinus1*nV[1]                               + nV[1]*dpsdUL[4] + rhouSL[1]*dSMdUL[4];
  
  AL[3][0] =    qL*uL[2]       - nV[2]*psiL*gammaMinus1/2.0   + nV[2]*dpsdUL[0] + rhouSL[2]*dSMdUL[0];
  AL[3][1] =     - uL[2]*nV[0] + nV[2]*gammaMinus1*uL[0]      + nV[2]*dpsdUL[1] + rhouSL[2]*dSMdUL[1];
  AL[3][2] =     - uL[2]*nV[1] + nV[2]*gammaMinus1*uL[1]      + nV[2]*dpsdUL[2] + rhouSL[2]*dSMdUL[2];
  AL[3][3] =  SL - qL          + nV[2]*(gamma-2.0)*uL[2]      + nV[2]*dpsdUL[3] + rhouSL[2]*dSMdUL[3];
  AL[3][4] = -gammaMinus1*nV[2]                               + nV[2]*dpsdUL[4] + rhouSL[2]*dSMdUL[4];
  
  AL[4][0] =      qL*(eL+pL)/Density_i - qL*psiL*(gamma-1.0)/2.0   + SM*dpsdUL[0] + (pS+eSL)*dSMdUL[0];
  AL[4][1] = - nV[0]*(eL+pL)/Density_i + gammaMinus1*uL[0]*qL      + SM*dpsdUL[1] + (pS+eSL)*dSMdUL[1];
  AL[4][2] = - nV[1]*(eL+pL)/Density_i + gammaMinus1*uL[1]*qL      + SM*dpsdUL[2] + (pS+eSL)*dSMdUL[2];
  AL[4][3] = - nV[2]*(eL+pL)/Density_i + gammaMinus1*uL[2]*qL      + SM*dpsdUL[3] + (pS+eSL)*dSMdUL[3];
  AL[4][4] =   SL-qL*gamma                                    + SM*dpsdUL[4] + (pS+eSL)*dSMdUL[4];
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      AL[i][j] *= omL;
  
  for (iVar = 0; iVar < nVar; iVar++)    AR[0][i] = omL*rhoSL*dSMdUR[i];
  for (iVar = 0; iVar < nVar; iVar++)    AR[1][i] = omL*(nV[0]*dpsdUR[i]+rhouSL[0]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[2][i] = omL*(nV[1]*dpsdUR[i]+rhouSL[1]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[3][i] = omL*(nV[2]*dpsdUR[i]+rhouSL[2]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[4][i] = omL*(dpsdUR[i]*SM+(pS+eSL)*dSMdUR[i]);
  
}

void UgpWithCvCompFlow::calcJacobianA(double (*A)[5], const double *vel, double pp, double rrho, const double *nV, double gamma, double surfVeloc) // nV is not normalized
{
 
  double kapm1 = (gamma - 1.0);
  
  double nVel[3];
  nVel[0] = vel[0]*nV[0];
  nVel[1] = vel[1]*nV[1];
  nVel[2] = vel[2]*nV[2];
  double U_k = nVel[0]+nVel[1]+nVel[2];
  double vSquHlf = 0.5*vecDotVec3d(vel, vel);
  double c = sqrt(gamma*pp/rrho);
  double inv_kap_m1 = 1.0/kapm1;
  
  A[0][0] =-surfVeloc;
  A[0][1] = nV[0];
  A[0][2] = nV[1];
  A[0][3] = nV[2];
  A[0][4] = 0.0;
  
  A[1][0] = -vel[0]*(nVel[1]+nVel[2])+nV[0]*(kapm1*vSquHlf-vel[0]*vel[0]);
  A[1][1] = (2.-gamma)*nVel[0]+U_k-surfVeloc;
  A[1][2] = vel[0]*nV[1]-kapm1*vel[1]*nV[0];
  A[1][3] = vel[0]*nV[2]-kapm1*vel[2]*nV[0];
  A[1][4] = kapm1*nV[0];
  
  A[2][0] = -vel[1]*(nVel[0]+nVel[2])+nV[1]*(kapm1*vSquHlf-vel[1]*vel[1]);
  A[2][1] = -kapm1*vel[0]*nV[1]+ vel[1]*nV[0];
  A[2][2] = (2.-gamma)*nVel[1]+U_k-surfVeloc;
  A[2][3] = vel[1]*nV[2]-kapm1*vel[2]*nV[1];
  A[2][4] = kapm1*nV[1];
  
  A[3][0] = -vel[2]*(nVel[0]+nVel[1])+nV[2]*(kapm1*vSquHlf-vel[2]*vel[2]);
  A[3][1] = -kapm1*vel[0]*nV[2]+vel[2]*nV[0];
  A[3][2] = -kapm1*vel[1]*nV[2]+vel[2]*nV[1];
  A[3][3] = (2.-gamma)*nVel[2]+U_k-surfVeloc;
  A[3][4] = kapm1*nV[2];
  
  A[4][0] = U_k*((gamma-2.)*vSquHlf-c*c*inv_kap_m1);
  A[4][1] = c*c*inv_kap_m1*nV[0]-kapm1*vel[0]*(nVel[1]+nVel[2])-(kapm1*vel[0]*vel[0]-vSquHlf)*nV[0];
  A[4][2] = c*c*inv_kap_m1*nV[1]-kapm1*vel[1]*(nVel[0]+nVel[2])-(kapm1*vel[1]*vel[1]-vSquHlf)*nV[1];
  A[4][3] = c*c*inv_kap_m1*nV[2]-kapm1*vel[2]*(nVel[0]+nVel[1])-(kapm1*vel[2]*vel[2]-vSquHlf)*nV[2];
  A[4][4] = gamma*U_k-surfVeloc;
  
}


#endif


CUpwRoe_Flow::CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  kappa = config->GetRoe_Kappa(); // 1 is unstable

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  delta_vel  = new double [nDim];
  delta_wave = new double [nVar];
  ProjFlux_i = new double [nVar];
  ProjFlux_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
}

CUpwRoe_Flow::~CUpwRoe_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwRoe_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};
  double ProjGridVel = 0.0;
  
  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));
  
  /*--- Primitive variables at point j ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));
  
  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  
  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  
  RoeSoundSpeed2 = (Gamma-1)*(RoeEnthalpy-0.5*sq_vel);
  
  /*--- Negative RoeSoundSpeed2, the jump 
   variables is too large, exit the subrotuine 
   without computing the fluxes ---*/
  
  if (RoeSoundSpeed2 <= 0.0) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][iVar] = 0.0;
        val_Jacobian_j[iVar][iVar] = 0.0;
      }
    }
    return;
  }

  RoeSoundSpeed = sqrt(RoeSoundSpeed2);
  
  /*--- Compute ProjFlux_i ---*/
  
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/
  
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);
  
  /*--- Compute P and Lambda (do it with the Normal) ---*/
  
  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Projected velocity adjustment due to mesh motion ---*/
  
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  /*--- Flow eigenvalues and entropy correctors ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
  /*--- Compute absolute value with Mavriplis' entropy correction ---*/
  
  MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;
  Delta = config->GetEntropyFix_Coeff();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Lambda[iVar] = max(fabs(Lambda[iVar]), Delta*MaxLambda);
  }
  
  if (!implicit) {
    
    /*--- Compute wave amplitudes (characteristics) ---*/
    
    proj_delta_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
      proj_delta_vel += delta_vel[iDim]*Normal[iDim];
    }
    delta_p = Pressure_j - Pressure_i;
    delta_rho = Density_j - Density_i;
    proj_delta_vel = proj_delta_vel/Area;
    
    if (nDim == 2) {
      delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
      delta_wave[1] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
      delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
      delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    } else {
      delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
      delta_wave[1] = UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0];
      delta_wave[2] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
      delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
      delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    }
    
    /*--- Roe's Flux approximation ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++)
        val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
    }
    
    /*--- Flux contribution due to grid motion ---*/
    
    if (grid_movement) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      }
    }
  }
  
  else {
    
    /*--- Compute inverse P ---*/
    
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobians of the inviscid flux, scaled by
     kappa because val_resconv ~ kappa*(fc_i+fc_j)*Normal ---*/
    
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, kappa, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, kappa, val_Jacobian_j);
    
    /*--- Diference variables iPoint and jPoint ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    
    /*--- Roe's Flux approximation ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      val_residual[iVar] = kappa*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
        
        val_residual[iVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
        val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        
      }
      
    }
    
    /*--- Jacobian contributions due to grid motion ---*/
    
    if (grid_movement) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
        
        /*--- Implicit terms ---*/
        
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
    
  }
  
}


CUpwGeneralRoe_Flow::CUpwGeneralRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();


  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  delta_vel  = new double [nDim];
  delta_wave = new double [nVar];
  ProjFlux_i = new double [nVar];
  ProjFlux_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
}

CUpwGeneralRoe_Flow::~CUpwGeneralRoe_Flow(void) {

  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;

}

void CUpwGeneralRoe_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

  double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

	/*--- Face area (norm or the normal vector) ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
	Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

	/*-- Unit Normal ---*/

  for (iDim = 0; iDim < nDim; iDim++)
		UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables at point i ---*/

  Velocity2_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
	  Velocity_i[iDim] = V_i[iDim+1];
	  Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
  }

  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  StaticEnthalpy_i = Enthalpy_i - 0.5*Velocity2_i;
  StaticEnergy_i = StaticEnthalpy_i - Pressure_i/Density_i;

  Kappa_i = S_i[1]/Density_i;
  Chi_i = S_i[0] - Kappa_i*StaticEnergy_i;
  SoundSpeed_i = sqrt(Chi_i + StaticEnthalpy_i*Kappa_i);

  /*--- Primitive variables at point j ---*/


  Velocity2_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    Velocity2_j += Velocity_j[iDim]*Velocity_j[iDim];

  }

  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;

  StaticEnthalpy_j = Enthalpy_j - 0.5*Velocity2_j;
  StaticEnergy_j = StaticEnthalpy_j - Pressure_j/Density_j;

  Kappa_j = S_j[1]/Density_j;
  Chi_j = S_j[0] - Kappa_j*StaticEnergy_j;
  SoundSpeed_j = sqrt(Chi_j + StaticEnthalpy_j*Kappa_j);

  /*--- Recompute conservative variables ---*/

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;

//	/*--- Roe-averaged variables at interface between i & j ---*/

    ComputeRoeAverage();


	/*--- Compute ProjFlux_i ---*/
	GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);

	/*--- Compute ProjFlux_j ---*/
	GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);

	/*--- Compute P and Lambda (do it with the Normal) ---*/

	GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, &RoeEnthalpy, &RoeChi, &RoeKappa, UnitNormal, P_Tensor);

	ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
		ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
		ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
	}

	/*--- Projected velocity adjustment due to mesh motion ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
		}
		ProjVelocity   -= ProjGridVel;
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j -= ProjGridVel;
	}

	/*--- Flow eigenvalues and entropy correctors ---*/
	for (iDim = 0; iDim < nDim; iDim++)
		Lambda[iDim] = ProjVelocity;

	Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
	Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;

//	/*--- Harten and Hyman (1983) entropy correction ---*/
//	for (iDim = 0; iDim < nDim; iDim++)
//		Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i, ProjVelocity_j-Lambda[iDim]));
//
//	Epsilon[nVar-2] = 4.0*max(0.0, max(Lambda[nVar-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar-2]));
//	Epsilon[nVar-1] = 4.0*max(0.0, max(Lambda[nVar-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar-1]));
//
//	for (iVar = 0; iVar < nVar; iVar++)
//		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
//			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
//		else
//			Lambda[iVar] = fabs(Lambda[iVar]);

	for (iVar = 0; iVar < nVar; iVar++)
		Lambda[iVar] = fabs(Lambda[iVar]);

	if (!implicit) {

		/*--- Compute wave amplitudes (characteristics) ---*/
		proj_delta_vel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
			proj_delta_vel += delta_vel[iDim]*Normal[iDim];
		}
		delta_p = Pressure_j - Pressure_i;
		delta_rho = Density_j - Density_i;
		proj_delta_vel = proj_delta_vel/Area;

		if (nDim == 2) {
			delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
			delta_wave[1] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
			delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
			delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
		} else {
			delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
			delta_wave[1] = UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0];
			delta_wave[2] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
			delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
			delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
		}

		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
			for (jVar = 0; jVar < nVar; jVar++)
				val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
		}

		/*--- Flux contribution due to grid motion ---*/
		if (grid_movement) {
			ProjVelocity = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
			}
		}
	}
	else {

		/*--- Compute inverse P ---*/

		GetPMatrix_inv(invP_Tensor, &RoeDensity, RoeVelocity, &RoeSoundSpeed, &RoeChi , &RoeKappa, UnitNormal);

		/*--- Jacobians of the inviscid flux, scaled by
        0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
		GetInviscidProjJac(Velocity_i, &Enthalpy_i, &Chi_i, &Kappa_i, Normal, 0.5, val_Jacobian_i);

		GetInviscidProjJac(Velocity_j, &Enthalpy_j, &Chi_j, &Kappa_j, Normal, 0.5, val_Jacobian_j);


		/*--- Diference variables iPoint and jPoint ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			Diff_U[iVar] = U_j[iVar]-U_i[iVar];

		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
			for (jVar = 0; jVar < nVar; jVar++) {
				Proj_ModJac_Tensor_ij = 0.0;
				/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
				for (kVar = 0; kVar < nVar; kVar++)
					Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
				val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
				val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
				val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
			}
		}

		/*--- Jacobian contributions due to grid motion ---*/
		if (grid_movement) {
			ProjVelocity = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
				/*--- Implicit terms ---*/
				val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
				val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
			}
		}

	}

}


void CUpwGeneralRoe_Flow::ComputeRoeAverage() {

	double delta_rhoStaticEnergy, err_P, s, D;//, stateSeparationLimit;
	// double tol = 10-6;
	//
	R = sqrt(fabs(Density_j/Density_i));
	RoeDensity = R*Density_i;
	sq_vel = 0;	for (iDim = 0; iDim < nDim; iDim++) {
		RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
		sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
	}

	RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
	delta_rho = Density_j - Density_i;
	delta_p = Pressure_j - Pressure_i;
	RoeKappa = 0.5*(Kappa_i + Kappa_j);
	RoeKappa = (Kappa_i + Kappa_j + 4*RoeKappa)/6;
	RoeChi = 0.5*(Chi_i + Chi_j);
	RoeChi = (Chi_i + Chi_j + 4*RoeChi)/6;

	//

	RoeKappaStaticEnthalpy = 0.5*(StaticEnthalpy_i*Kappa_i + StaticEnthalpy_j*Kappa_j);
	RoeKappaStaticEnthalpy = (StaticEnthalpy_i*Kappa_i + StaticEnthalpy_j*Kappa_j + 4*RoeKappaStaticEnthalpy)/6;
	s = RoeChi + RoeKappaStaticEnthalpy;
	D = s*s*delta_rho*delta_rho + delta_p*delta_p;
	delta_rhoStaticEnergy = Density_j*StaticEnergy_j - Density_i*StaticEnergy_i;
	err_P = delta_p - RoeChi*delta_rho - RoeKappa*delta_rhoStaticEnergy;


	if (abs((D - delta_p*err_P)/Density_i)>1e-3 && abs(delta_rho/Density_i)>1e-3 && s/Density_i > 1e-3) {

		RoeKappa = (D*RoeKappa)/(D - delta_p*err_P);
		RoeChi = (D*RoeChi+ s*s*delta_rho*err_P)/(D - delta_p*err_P);

	}

	RoeSoundSpeed = sqrt(RoeChi + RoeKappa*(RoeEnthalpy-0.5*sq_vel));
}

CUpwMSW_Flow::CUpwMSW_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  /*--- Set booleans from CConfig settings ---*/
	implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  
  /*--- Allocate arrays ---*/
	Diff_U   = new double [nVar];
  Fc_i	   = new double [nVar];
	Fc_j	   = new double [nVar];
	Lambda_i = new double [nVar];
  Lambda_j = new double [nVar];
  
	u_i		   = new double [nDim];
	u_j		   = new double [nDim];
  ust_i    = new double [nDim];
  ust_j    = new double [nDim];
  Vst_i    = new double [nPrimVar];
  Vst_j    = new double [nPrimVar];
  Ust_i    = new double [nVar];
  Ust_j    = new double [nVar];
  
  Velst_i    = new double [nDim];
  Velst_j    = new double [nDim];
  
	P_Tensor		= new double* [nVar];
	invP_Tensor	= new double* [nVar];
	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar]    = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
	}
  
}

CUpwMSW_Flow::~CUpwMSW_Flow(void) {
  
	delete [] Diff_U;
  delete [] Fc_i;
	delete [] Fc_j;
	delete [] Lambda_i;
  delete [] Lambda_j;
  
  delete [] u_i;
  delete [] u_j;
  delete [] ust_i;
  delete [] ust_j;
  delete [] Ust_i;
  delete [] Vst_i;
  delete [] Ust_j;
  delete [] Vst_j;
  delete [] Velst_i;
  delete [] Velst_j;
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;

}

void CUpwMSW_Flow::ComputeResidual(double *val_residual,
                                   double **val_Jacobian_i,
                                   double **val_Jacobian_j, CConfig *config) {
  
	unsigned short iDim, iVar, jVar, kVar;
  double P_i, P_j;
  double ProjVel_i, ProjVel_j, ProjVelst_i, ProjVelst_j;
  double sqvel_i, sqvel_j;
	double alpha, w, dp, onemw;
  double Proj_ModJac_Tensor_i, Proj_ModJac_Tensor_j;
  
  /*--- Set parameters in the numerical method ---*/
  alpha = 6.0;
  
  /*--- Calculate supporting geometry parameters ---*/
  
	Area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Initialize flux & Jacobian vectors ---*/
  
	for (iVar = 0; iVar < nVar; iVar++) {
		Fc_i[iVar] = 0.0;
		Fc_j[iVar] = 0.0;
	}
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][jVar] = 0.0;
        val_Jacobian_j[iVar][jVar] = 0.0;
      }
    }
  }
  
  /*--- Load variables from nodes i & j ---*/
  
  rhos_i = V_i[0];
  rhos_j = V_j[0];
  for (iDim = 0; iDim < nDim; iDim++) {
    u_i[iDim] = V_i[iDim+1];
    u_j[iDim] = V_j[iDim+1];
  }
  P_i = V_i[nDim+1];
  P_j = V_j[nDim+1];
  
  /*--- Calculate supporting quantities ---*/
  
  sqvel_i   = 0.0; sqvel_j   = 0.0;
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel_i   += u_i[iDim]*u_i[iDim];
    sqvel_j   += u_j[iDim]*u_j[iDim];
    ProjVel_i += u_i[iDim]*UnitNormal[iDim];
    ProjVel_j += u_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Calculate the state weighting function ---*/
  
  dp = fabs(P_j-P_i) / min(P_j, P_i);
  w = 0.5 * (1.0/(pow(alpha*dp,2.0) +1.0));
  onemw = 1.0 - w;
  
  /*--- Calculate weighted state vector (*) for i & j ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Ust_i[iVar] = onemw*U_i[iVar] + w*U_j[iVar];
    Ust_j[iVar] = onemw*U_j[iVar] + w*U_i[iVar];
  }
  for (iVar = 0; iVar < nDim+5; iVar++) {
    Vst_i[iVar] = onemw*V_i[iVar] + w*V_j[iVar];
    Vst_j[iVar] = onemw*V_j[iVar] + w*V_i[iVar];
  }
  ProjVelst_i = onemw*ProjVel_i + w*ProjVel_j;
  ProjVelst_j = onemw*ProjVel_j + w*ProjVel_i;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Velst_i[iDim] = Vst_i[iDim+1];
    Velst_j[iDim] = Vst_j[iDim+1];
  }
  
  /*--- Flow eigenvalues at i (Lambda+) --- */
  
  for (iDim = 0; iDim < nDim; iDim++) {
  Lambda_i[iDim]      = 0.5*(ProjVelst_i + fabs(ProjVelst_i));
  }

  Lambda_i[nDim] = 0.5*( ProjVelst_i + Vst_i[nDim+4] + fabs(ProjVelst_i + Vst_i[nDim+4])  );
  Lambda_i[nDim+1]   = 0.5*( ProjVelst_i - Vst_i[nDim+4] + fabs(ProjVelst_i - Vst_i[nDim+4])  );
  
  /*--- Compute projected P, invP, and Lambda ---*/
  
  GetPMatrix(&Vst_i[nDim+2], Velst_i, &Vst_i[nDim+4], UnitNormal, P_Tensor);
  GetPMatrix_inv(&Vst_i[nDim+2], Velst_i, &Vst_i[nDim+4], UnitNormal, invP_Tensor);
  
  /*--- Projected flux (f+) at i ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_i = 0.0;
      
      /*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
      
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_i += P_Tensor[iVar][kVar]*Lambda_i[kVar]*invP_Tensor[kVar][jVar];
      Fc_i[iVar] += Proj_ModJac_Tensor_i*U_i[jVar]*Area;
      if (implicit)
        val_Jacobian_i[iVar][jVar] += Proj_ModJac_Tensor_i*Area;
    }
  }
  
	/*--- Flow eigenvalues at j (Lambda-) ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Lambda_j[iDim]          = 0.5*(ProjVelst_j - fabs(ProjVelst_j));
  }
  Lambda_j[nDim] = 0.5*(     ProjVelst_j + Vst_j[nDim+4] -
                                   fabs(ProjVelst_j + Vst_j[nDim+4])  );
  Lambda_j[nDim+1]   = 0.5*(     ProjVelst_j - Vst_j[nDim+4] -
                                   fabs(ProjVelst_j - Vst_j[nDim+4])  );
  
  /*--- Compute projected P, invP, and Lambda ---*/
  
  GetPMatrix(&Vst_j[nDim+2], Velst_j, &Vst_j[nDim+4], UnitNormal, P_Tensor);
  GetPMatrix_inv(&Vst_j[nDim+2], Velst_j, &Vst_j[nDim+4], UnitNormal, invP_Tensor);
  
	/*--- Projected flux (f-) ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_j = 0.0;
      /*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_j += P_Tensor[iVar][kVar]*Lambda_j[kVar]*invP_Tensor[kVar][jVar];
      Fc_j[iVar] += Proj_ModJac_Tensor_j*U_j[jVar]*Area;
      if (implicit)
        val_Jacobian_j[iVar][jVar] += Proj_ModJac_Tensor_j*Area;
    }
  }
  
	/*--- Flux splitting ---*/
  
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual[iVar] = Fc_i[iVar]+Fc_j[iVar];
	}
  
}

CUpwTurkel_Flow::CUpwTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Beta_min = config->GetminTurkelBeta();
  Beta_max = config->GetmaxTurkelBeta();
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  ProjFlux_i = new double [nVar];
  ProjFlux_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  absPeJac = new double* [nVar];
  invRinvPe = new double* [nVar];
  R_Tensor  = new double* [nVar];
  Matrix    = new double* [nVar];
  Art_Visc  = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    absPeJac[iVar] = new double [nVar];
    invRinvPe[iVar] = new double [nVar];
    Matrix[iVar] = new double [nVar];
    Art_Visc[iVar] = new double [nVar];
    R_Tensor[iVar] = new double [nVar];
  }
}

CUpwTurkel_Flow::~CUpwTurkel_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] absPeJac[iVar];
    delete [] invRinvPe[iVar];
    delete [] Matrix[iVar];
    delete [] Art_Visc[iVar];
    delete [] R_Tensor[iVar];
  }
  delete [] Matrix;
  delete [] Art_Visc;
  delete [] absPeJac;
  delete [] invRinvPe;
  delete [] R_Tensor;
  
}

void CUpwTurkel_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));

  /*--- Primitive variables at point j ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));

  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  
  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
  RoePressure = RoeDensity/Gamma*RoeSoundSpeed*RoeSoundSpeed;
  
  /*--- Compute ProjFlux_i ---*/
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Projected velocity adjustment due to mesh motion ---*/
  if (grid_movement) {
    double ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  /*--- First few flow eigenvalues of A.Normal with the normal---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  local_Mach = sqrt(sq_vel)/RoeSoundSpeed;
  Beta 	   = max(Beta_min, min(local_Mach, Beta_max));
  Beta2 	   = Beta*Beta;
  
  one_m_Betasqr 		   = 1.0 - Beta2;  // 1-Beta*Beta
  one_p_Betasqr 		   = 1.0 + Beta2;  // 1+Beta*Beta
  sqr_one_m_Betasqr_Lam1 = pow((one_m_Betasqr*Lambda[0]),2); // [(1-Beta^2)*Lambda[0]]^2
  sqr_two_Beta_c_Area    = pow(2.0*Beta*RoeSoundSpeed*Area,2); // [2*Beta*c*Area]^2
  
  /*--- The rest of the flow eigenvalues of preconditioned matrix---*/
  Lambda[nVar-2] = 0.5 * ( one_p_Betasqr*Lambda[0] + sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
  Lambda[nVar-1] = 0.5 * ( one_p_Betasqr*Lambda[0] - sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
  
  s_hat = 1.0/Area * (Lambda[nVar-1] - Lambda[0]*Beta2);
  r_hat = 1.0/Area * (Lambda[nVar-2] - Lambda[0]*Beta2);
  t_hat = 0.5/Area * (Lambda[nVar-1] - Lambda[nVar-2]);
  rhoB2a2 = RoeDensity*Beta2*RoeSoundSpeed*RoeSoundSpeed;
  
  /*--- Diference variables iPoint and jPoint and absolute value of the eigen values---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    Lambda[iVar] = fabs(Lambda[iVar]);
  }
  
  /*--- Compute the absolute Preconditioned Jacobian in entropic Variables (do it with the Unitary Normal) ---*/
  GetPrecondJacobian(Beta2, r_hat, s_hat, t_hat, rhoB2a2, Lambda, UnitNormal, absPeJac);
  
  /*--- Compute the matrix from entropic variables to conserved variables ---*/
  GetinvRinvPe(Beta2, RoeEnthalpy, RoeSoundSpeed, RoeDensity, RoeVelocity, invRinvPe);
  
  /*--- Compute the matrix from entropic variables to conserved variables ---*/
  GetRMatrix(RoePressure, RoeSoundSpeed, RoeDensity, RoeVelocity, R_Tensor);
  
  if (implicit) {
    /*--- Jacobians of the inviscid flux, scaled by
     0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
  }
  
  for (iVar = 0; iVar < nVar; iVar ++) {
    for (jVar = 0; jVar < nVar; jVar ++) {
      Matrix[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Matrix[iVar][jVar]  += absPeJac[iVar][kVar]*R_Tensor[kVar][jVar];
    }
  }
  
  for (iVar = 0; iVar < nVar; iVar ++) {
    for (jVar = 0; jVar < nVar; jVar ++) {
      Art_Visc[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Art_Visc[iVar][jVar]  += invRinvPe[iVar][kVar]*Matrix[kVar][jVar];
    }
  }
  
  /*--- Roe's Flux approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual[iVar] -= 0.5*Art_Visc[iVar][jVar]*Diff_U[jVar];
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += 0.5*Art_Visc[iVar][jVar];
        val_Jacobian_j[iVar][jVar] -= 0.5*Art_Visc[iVar][jVar];
      }
    }
  }
  
  /*--- Contributions due to mesh motion---*/
  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      /*--- Implicit terms ---*/
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }
  
}


CUpwArtComp_Flow::CUpwArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  gravity = config->GetGravityForce();
  Froude = config->GetFroude();
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  ProjFlux_i = new double [nVar];
  ProjFlux_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
  
}

CUpwArtComp_Flow::~CUpwArtComp_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwArtComp_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*--- Compute and unitary normal vector ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = Normal[iDim]/Area;
    if (fabs(UnitNormal[iDim]) < EPS) UnitNormal[iDim] = EPS;
  }
  
  /*--- Set velocity and pressure variables at points iPoint and jPoint ---*/
  
  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+1];  DensityInc_j = V_j[nDim+1];
  BetaInc2_i =    V_i[nDim+2];  BetaInc2_j = V_j[nDim+2];
  
  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
    ProjVelocity += MeanVelocity[iDim]*Normal[iDim];
  }
  
  /*--- Mean variables at points iPoint and jPoint ---*/
  
  MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
  MeanPressure = 0.5*(Pressure_i + Pressure_j);
  MeanBetaInc2 = 0.5*(BetaInc2_i + BetaInc2_j);
  MeanSoundSpeed = sqrt(ProjVelocity*ProjVelocity + (MeanBetaInc2/MeanDensity) * Area * Area);
  
  /*--- Compute ProjFlux_i ---*/
  
  GetInviscidArtCompProjFlux(&DensityInc_i, Velocity_i, &Pressure_i, &BetaInc2_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/
  
  GetInviscidArtCompProjFlux(&DensityInc_j, Velocity_j, &Pressure_j, &BetaInc2_j, Normal, ProjFlux_j);
  
  /*--- Compute P and Lambda (matrix of eigenvalues) ---*/
  
  GetPArtCompMatrix(&MeanDensity, MeanVelocity, &MeanBetaInc2, UnitNormal, P_Tensor);
  
  /*--- Flow eigenvalues ---*/
  
  if (nDim == 2) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity + MeanSoundSpeed;
    Lambda[2] = ProjVelocity - MeanSoundSpeed;
  }
  if (nDim == 3) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity;
    Lambda[2] = ProjVelocity + MeanSoundSpeed;
    Lambda[3] = ProjVelocity - MeanSoundSpeed;
  }
  
  /*--- Absolute value of the eigenvalues ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  /*--- Compute inverse P ---*/
  
  GetPArtCompMatrix_inv(&MeanDensity, MeanVelocity, &MeanBetaInc2, UnitNormal, invP_Tensor);

  /*--- Jacobian of the inviscid flux ---*/

  if (implicit) {
    GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, Normal, 0.5, val_Jacobian_j);
  }
  
  /*--- Diference variables iPoint and jPoint ---*/
  
  Diff_U[0] = Pressure_j - Pressure_i;
  for (iDim = 0; iDim < nDim; iDim++)
    Diff_U[iDim+1] = Velocity_j[iDim]*DensityInc_i - Velocity_i[iDim]*DensityInc_j;
  
  /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
      val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar];
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij;
      }
    }
  }
  
}

CUpwArtComp_FreeSurf_Flow::CUpwArtComp_FreeSurf_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  gravity = config->GetGravityForce();
  Froude = config->GetFroude();
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  ProjFlux_i = new double [nVar];
  ProjFlux_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
  
}

CUpwArtComp_FreeSurf_Flow::~CUpwArtComp_FreeSurf_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwArtComp_FreeSurf_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Compute face area ---*/
  
  Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*--- Compute and unitary normal vector ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = Normal[iDim]/Area;
    if (fabs(UnitNormal[iDim]) < EPS) UnitNormal[iDim] = EPS;
  }
  
  /*--- Set velocity and pressure and level set variables at points iPoint and jPoint ---*/
  
  Pressure_i = V_i[0]; Pressure_j = V_j[0];
  DensityInc_i = V_i[nDim+1]; DensityInc_j = V_j[nDim+1];
  BetaInc2_i = V_i[nDim+2]; BetaInc2_j = V_j[nDim+2];
  LevelSet_i = V_i[nDim+5]; LevelSet_j = V_j[nDim+5];
  Distance_i = V_i[nDim+6]; Distance_j = V_j[nDim+6];

  if (fabs(LevelSet_i) < EPS) LevelSet_i = EPS;
  if (fabs(LevelSet_j) < EPS) LevelSet_j = EPS;
  
  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1]; if (fabs(Velocity_i[iDim]) < EPS) Velocity_i[iDim] = EPS;
    Velocity_j[iDim] = V_j[iDim+1]; if (fabs(Velocity_j[iDim]) < EPS) Velocity_j[iDim] = EPS;
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
    ProjVelocity += MeanVelocity[iDim]*Normal[iDim];
  }
  
  double epsilon = config->GetFreeSurface_Thickness(), Delta = 0.0;
  if (LevelSet_i < -epsilon) Delta = 0.0;
  if (fabs(LevelSet_i) <= epsilon) Delta = 0.5*(1.0+cos(PI_NUMBER*LevelSet_i/epsilon))/epsilon;
  if (LevelSet_i > epsilon) Delta = 0.0;
  dDensityInc_i = (1.0 - config->GetRatioDensity())*Delta*config->GetDensity_FreeStreamND();
  
  if (LevelSet_j < -epsilon) Delta = 0.0;
  if (fabs(LevelSet_j) <= epsilon) Delta = 0.5*(1.0+cos(PI_NUMBER*LevelSet_j/epsilon))/epsilon;
  if (LevelSet_j > epsilon) Delta = 0.0;
  dDensityInc_j = (1.0 - config->GetRatioDensity())*Delta*config->GetDensity_FreeStreamND();
  
  /*--- Mean variables at points iPoint and jPoint ---*/
  
  MeanDensityInc   = 0.5*(DensityInc_i + DensityInc_j);
  dMeanDensityInc   = 0.5*(dDensityInc_i + dDensityInc_j);
  MeanPressure  = 0.5*(Pressure_i + Pressure_j);
  MeanLevelSet  = 0.5*(LevelSet_i + LevelSet_j);
  MeanBetaInc2  = 0.5*(BetaInc2_i + BetaInc2_j);
  
  /*--- Compute ProjFlux_i ---*/
  
  GetInviscidArtComp_FreeSurf_ProjFlux(&DensityInc_i, Velocity_i, &Pressure_i, &BetaInc2_i, &LevelSet_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/
  
  GetInviscidArtComp_FreeSurf_ProjFlux(&DensityInc_j, Velocity_j, &Pressure_j, &BetaInc2_j, &LevelSet_j, Normal, ProjFlux_j);
  
  /*--- Compute P and Lambda (matrix of eigenvalues) ---*/
  
  GetPArtComp_FreeSurf_Matrix(&MeanDensityInc, &dMeanDensityInc, MeanVelocity, &MeanBetaInc2, &MeanLevelSet, UnitNormal, P_Tensor);
  
  /*--- Flow eigenvalues ---*/
  
  double a = MeanBetaInc2/MeanDensityInc, b = MeanLevelSet/MeanDensityInc, c = dMeanDensityInc;
  double e = (2.0 - b*c)*ProjVelocity, f = sqrt(4.0*a*Area*Area + e*e);
  
  if (nDim == 2) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity;
    Lambda[2] = 0.5*(e - f);
    Lambda[3] = 0.5*(e + f);
  }
  if (nDim == 3) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity;
    Lambda[2] = ProjVelocity;
    Lambda[3] = 0.5*(e - f);
    Lambda[4] = 0.5*(e + f);
  }
  
  /*--- Absolute value of the eigenvalues ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  /*--- Compute inverse P ---*/
  
  GetPArtComp_FreeSurf_Matrix_inv(&MeanDensityInc, &dMeanDensityInc, MeanVelocity, &MeanBetaInc2, &MeanLevelSet, UnitNormal, invP_Tensor);

  /*--- Jacobian of the inviscid flux ---*/

  if (implicit) {
    GetInviscidArtComp_FreeSurf_ProjJac(&DensityInc_i, &dDensityInc_i, Velocity_i, &BetaInc2_i, &LevelSet_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidArtComp_FreeSurf_ProjJac(&DensityInc_j, &dDensityInc_j, Velocity_j, &BetaInc2_j, &LevelSet_j, Normal, 0.5, val_Jacobian_j);
  }
  
  /*--- Diference of conservative iPoint and jPoint ---*/
  
  Diff_U[0] = V_j[0] - V_i[0];
  for (iDim = 0; iDim < nDim; iDim++)
    Diff_U[iDim+1] = Velocity_j[iDim]*DensityInc_j - Velocity_i[iDim]*DensityInc_i;
  Diff_U[nDim] = LevelSet_j - LevelSet_i;

  /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      
      Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
      val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar];
      
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij;
      }
      
    }
  }
  
}

CCentJSTArtComp_Flow::CCentJSTArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  grid_movement = config->GetGrid_Movement();
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  gravity = config->GetGravityForce();
  Froude = config->GetFroude();
  
  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new double [nVar];
  Diff_Lapl = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  ProjFlux = new double [nVar];
  
}

CCentJSTArtComp_Flow::~CCentJSTArtComp_Flow(void) {
  
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
  
}

void CCentJSTArtComp_Flow::ComputeResidual(double *val_residual,
                                           double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  double U_i[4] = {0.0,0.0,0.0,0.0}, U_j[4] = {0.0,0.0,0.0,0.0};

  /*--- Primitive variables at point i and j ---*/
  
  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+1];  DensityInc_j = V_j[nDim+1];
  BetaInc2_i =    V_i[nDim+2];  BetaInc2_j = V_j[nDim+2];

  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
  }
  
  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Pressure_i; U_j[0] = Pressure_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = DensityInc_i*Velocity_i[iDim]; U_j[iDim+1] = DensityInc_j*Velocity_j[iDim];
  }
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
  MeanPressure = 0.5*(Pressure_i + Pressure_j);
  MeanBetaInc2 = 0.5*(BetaInc2_i + BetaInc2_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidArtCompProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, Normal, ProjFlux);
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux ---*/
  
  if (implicit) {
    GetInviscidArtCompProjJac(&MeanDensity, MeanVelocity, &MeanBetaInc2, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Computes differences between Laplacians and conservative variables ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  }
  
  /*--- Compute the local espectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
  SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);
  
  Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
  Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
  sc4 = sc2*sc2/4.0;
  
  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
  
  if (implicit) {
    cte_0 = (Epsilon_2 + Epsilon_4*double(Neighbor_i+1))*StretchingFactor*MeanLambda;
    cte_1 = (Epsilon_2 + Epsilon_4*double(Neighbor_j+1))*StretchingFactor*MeanLambda;
    
    for (iVar = 0; iVar < nVar; iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_1;
    }
  }
  
}

CCentLaxArtComp_Flow::CCentLaxArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  gravity = config->GetGravityForce();
  Froude = config->GetFroude();
  
  /*--- Artificial dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  ProjFlux = new double [nVar];
  
}

CCentLaxArtComp_Flow::~CCentLaxArtComp_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
  
}

void CCentLaxArtComp_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                           CConfig *config) {
  
  double U_i[4] = {0.0,0.0,0.0,0.0}, U_j[4] = {0.0,0.0,0.0,0.0};

  /*--- Conservative variables at point i and j ---*/
  
  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+1];  DensityInc_j = V_j[nDim+1];
  BetaInc2_i =    V_i[nDim+2];  BetaInc2_j = V_j[nDim+2];
  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  
  /*--- Recompute conservative variables ---*/

  U_i[0] = Pressure_i; U_j[0] = Pressure_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = DensityInc_i*Velocity_i[iDim]; U_j[iDim+1] = DensityInc_j*Velocity_j[iDim];
  }
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(DensityInc_i+DensityInc_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanBetaInc2 = 0.5*(BetaInc2_i+BetaInc2_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidArtCompProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, Normal, ProjFlux);
  
  /*--- Compute inviscid residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux ---*/
  
  if (implicit) {
    GetInviscidArtCompProjJac(&MeanDensity, MeanVelocity, &MeanBetaInc2, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Computes differences btw. conservative variables ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  
  /*--- Compute the local espectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0; ProjVelocity_j = 0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
  SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);
  
  Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
  Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;
  MeanLambda = 0.5*(Local_Lambda_i + Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc0 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*double(nDim)/3.0;
  
  /*--- Compute viscous part of the residual ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_Jacobian_i[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
      val_Jacobian_j[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
    }
  }
  
}

CAvgGrad_Flow::CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  PrimVar_i = new double [nDim+3];
  PrimVar_j = new double [nDim+3];
  Mean_PrimVar = new double [nDim+3];
  
  Mean_GradPrimVar = new double* [nDim+1];
  for (iVar = 0; iVar < nDim+1; iVar++)
    Mean_GradPrimVar[iVar] = new double [nDim];
  
}

CAvgGrad_Flow::~CAvgGrad_Flow(void) {

  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  for (iVar = 0; iVar < nDim+1; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
}

void CAvgGrad_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  for (iVar = 0; iVar < nDim+3; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

	  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5]; Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6]; Eddy_Viscosity_j = V_j[nDim+6];

  /*--- Mean Viscosities and turbulent kinetic energy---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke = 0.5*(turb_ke_i + turb_ke_j);
  
  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = .5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
    }
  }

  /*--- Get projected flux tensor ---*/
  
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  
  if (implicit) {
    
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);
    
    if (dist_ij == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
	for (jVar = 0; jVar < nVar; jVar++) {
	  val_Jacobian_i[iVar][jVar] = 0.0;
	  val_Jacobian_j[iVar][jVar] = 0.0;
	}
      }
    }
    else {
      GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
			 dist_ij, UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}

CAvgGrad_Flow2::CAvgGrad_Flow2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  PrimVar_i = new double [nDim+3];
  PrimVar_j = new double [nDim+3];
  Mean_PrimVar = new double [nDim+3];
  
  Mean_GradPrimVar = new double* [nDim+1];
  for (iVar = 0; iVar < nDim+1; iVar++)
    Mean_GradPrimVar[iVar] = new double [nDim];

  Mean_Tau = new double* [nDim];
  for (iVar = 0; iVar < nDim; iVar++)
    Mean_Tau[iVar] = new double [nDim];
  
}

CAvgGrad_Flow2::~CAvgGrad_Flow2(void) {

  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  for (iVar = 0; iVar < nDim+1; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  for (iVar = 0; iVar < nDim; iVar++)
    delete [] Mean_Tau[iVar];
  delete [] Mean_Tau;
  
}

void CAvgGrad_Flow2::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  for (iVar = 0; iVar < nDim+3; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

	  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5]; Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6]; Eddy_Viscosity_j = V_j[nDim+6];

  /*--- Mean Viscosities and turbulent kinetic energy---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke = 0.5*(turb_ke_i + turb_ke_j);
  
  /*--- Mean gradient approximation ---*/
  iVar = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Mean_GradPrimVar[iVar][iDim] = .5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
  }
  
  /*--- Shear stress ---*/
  Total_Viscosity_i = Laminar_Viscosity_i+Eddy_Viscosity_i;
  Total_Viscosity_j = Laminar_Viscosity_j+Eddy_Viscosity_j;

  visc_div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    visc_div_vel += .5*(Total_Viscosity_i*PrimVar_Grad_i[iDim+1][iDim] + 
			Total_Viscosity_j*PrimVar_Grad_j[iDim+1][iDim]);
  
  for (iVar = 0; iVar < nDim; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_Tau[iVar][iDim] = 0.5*(Total_Viscosity_i*PrimVar_Grad_i[iVar+1][iDim] + 
				  Total_Viscosity_j*PrimVar_Grad_j[iVar+1][iDim] +
				  Total_Viscosity_i*PrimVar_Grad_i[iDim+1][iVar] + 
				  Total_Viscosity_j*PrimVar_Grad_j[iDim+1][iVar]
				  );
    }
  }

  for (iDim = 0 ; iDim < nDim; iDim++)
    Mean_Tau[iDim][iDim] -= TWO3*(visc_div_vel + .5*Mean_turb_ke);

  /*--- Get projected flux tensor ---*/
  
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_Tau, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  
  if (implicit) {
    
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);
    
    if (dist_ij == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
	for (jVar = 0; jVar < nVar; jVar++) {
	  val_Jacobian_i[iVar][jVar] = 0.0;
	  val_Jacobian_j[iVar][jVar] = 0.0;
	}
      }
    }
    else {
      GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
			 dist_ij, UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}

CGeneralAvgGrad_Flow::CGeneralAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Compressible flow, primitive variables nDim+3, (vx, vy, vz, P, rho, h) ---*/
  PrimVar_i = new double [nDim+4];
  PrimVar_j = new double [nDim+4];
  Mean_PrimVar = new double [nDim+4];
  Mean_SecVar = new double [2];
  
  /*--- Compressible flow, primitive gradient variables nDim+3, (T, vx, vy, vz) ---*/
  Mean_GradPrimVar = new double* [nDim+1];
  for (iVar = 0; iVar < nDim+1; iVar++)
    Mean_GradPrimVar[iVar] = new double [nDim];
}

CGeneralAvgGrad_Flow::~CGeneralAvgGrad_Flow(void) {
  
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  
  delete [] Mean_SecVar;

  for (iVar = 0; iVar < nDim+1; iVar++)
  delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
}

void CGeneralAvgGrad_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mean primitive variables ---*/
  for (iVar = 0; iVar < nDim+4; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  Laminar_Viscosity_i = V_i[nDim+5];    Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6];       Eddy_Viscosity_j = V_j[nDim+6];
  Thermal_Conductivity_i = V_i[nDim+7]; Thermal_Conductivity_j = V_j[nDim+7];
  Cp_i = V_i[nDim+8]; Cp_j = V_j[nDim+8];

  /*--- Mean secondary variables ---*/
  for (iVar = 0; iVar < 2; iVar++) {
    Mean_SecVar[iVar] = 0.5*(S_i[iVar+2]+S_j[iVar+2]);
  }
  
  /*--- Mean Viscosities and turbulent kinetic energy---*/
  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke              = 0.5*(turb_ke_i + turb_ke_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);
  Mean_Cp                   = 0.5*(Cp_i + Cp_j);

  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
    }
  }
  
  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux( Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
		              Mean_Thermal_Conductivity, Mean_Cp );
  
  /*--- Update viscous residual ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  if (implicit) {
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);
    
    if (dist_ij == 0.0) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
//        GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
//                           dist_ij, UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
        GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar, Mean_SecVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_Thermal_Conductivity, Mean_Cp,
                           dist_ij, UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}

CAvgGradArtComp_Flow::CAvgGradArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Incompressible flow, primitive variables nDim+1, (P, vx, vy, vz) ---*/
  
  Mean_GradPrimVar = new double* [nVar];
  
  /*--- Incompressible flow, gradient primitive variables nDim+1, (P, vx, vy, vz) ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPrimVar[iVar] = new double [nDim];
  
}

CAvgGradArtComp_Flow::~CAvgGradArtComp_Flow(void) {
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
}

void CAvgGradArtComp_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i,
                                           double **val_Jacobian_j, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
  Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  
  /*--- Mean Viscosities ---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  
  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
  
  /*--- Get projected flux tensor ---*/
  
  GetViscousArtCompProjFlux(Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Implicit part ---*/
  
  if (implicit) {
    
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);
    
    if (dist_ij == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal,
                                Area, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}

CAvgGradCorrected_Flow::CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  limiter = config->GetViscous_Limiter_Flow();

	PrimVar_i = new double [nDim+3];
	PrimVar_j = new double [nDim+3];
	Mean_PrimVar = new double [nDim+3];
  
	Proj_Mean_GradPrimVar_Edge = new double [nDim+1];
	Mean_GradPrimVar = new double* [nDim+1];
	for (iVar = 0; iVar < nDim+1; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];
	Edge_Vector = new double [nDim];
  
}
CAvgGradCorrected_Flow::~CAvgGradCorrected_Flow(void) {

	delete [] PrimVar_i;
	delete [] PrimVar_j;
	delete [] Mean_PrimVar;
	delete [] Proj_Mean_GradPrimVar_Edge;
	delete [] Edge_Vector;
	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
  
}
void CAvgGradCorrected_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	/*--- Normalized normal vector ---*/
  
	Area = 0.0;
	for (iDim = 0; iDim < nDim; iDim++)
		Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
  
	for (iDim = 0; iDim < nDim; iDim++)
		UnitNormal[iDim] = Normal[iDim]/Area;
  
  for (iVar = 0; iVar < nDim+3; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5]; Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6]; Eddy_Viscosity_j = V_j[nDim+6];
  
  /*--- Mean Viscosities and turbulent kinetic energy ---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke = 0.5*(turb_ke_i + turb_ke_j);
  
  /*--- Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nDim+1; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    
    if (!limiter) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
        Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
      }
    }
    else {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim]*PrimVar_Lim_i[iVar] +
					    PrimVar_Grad_j[iVar][iDim]*PrimVar_Lim_j[iVar]);
        Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
      }
    }
    

    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                         (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
      }
    }
  }

  /*--- Get projected flux tensor ---*/
  
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
  /*--- Save residual value ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  
  if (implicit) {
    
    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
	for (jVar = 0; jVar < nVar; jVar++) {
	  val_Jacobian_i[iVar][jVar] = 0.0;
	  val_Jacobian_j[iVar][jVar] = 0.0;
	}
      }
    }
    else {
      GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
			 sqrt(dist_ij_2), UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}

CAvgGradCorrected_Flow2::CAvgGradCorrected_Flow2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  limiter = config->GetViscous_Limiter_Flow();

	PrimVar_i = new double [nDim+3];
	PrimVar_j = new double [nDim+3];
	Mean_PrimVar = new double [nDim+3];
  
	Proj_Mean_GradPrimVar_Edge = new double [nDim+1];
	Mean_GradPrimVar = new double* [nDim+1];
	for (iVar = 0; iVar < nDim+1; iVar++)
		Mean_GradPrimVar[iVar] = new double [nDim];

	Mean_Tau = new double* [nDim];
	for (iVar = 0; iVar < nDim; iVar++)
		Mean_Tau[iVar] = new double [nDim];

	Edge_Vector = new double [nDim];
  
}
CAvgGradCorrected_Flow2::~CAvgGradCorrected_Flow2(void) {

	delete [] PrimVar_i;
	delete [] PrimVar_j;
	delete [] Mean_PrimVar;
	delete [] Proj_Mean_GradPrimVar_Edge;
	delete [] Edge_Vector;
	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] Mean_GradPrimVar[iVar];
	delete [] Mean_GradPrimVar;
	for (iVar = 0; iVar < nDim+1; iVar++)
		delete [] Mean_Tau[iVar];
	delete [] Mean_Tau;
  
}
void CAvgGradCorrected_Flow2::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  for (iVar = 0; iVar < nDim+3; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5]; Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6]; Eddy_Viscosity_j = V_j[nDim+6];
  
  /*--- Mean Viscosities and turbulent kinetic energy ---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke = 0.5*(turb_ke_i + turb_ke_j);
  
  /*--- Projection of the mean gradient in the direction of the edge ---*/
  iVar = 0;
  Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    
  if (!limiter) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim]*PrimVar_Lim_i[iVar] +
					  PrimVar_Grad_j[iVar][iDim]*PrimVar_Lim_j[iVar]);
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    }
  }
    
  if (dist_ij_2 != 0.0) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
				       (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
    }
  }

  /*--- Shear stress ---*/
  Total_Viscosity_i = Laminar_Viscosity_i+Eddy_Viscosity_i;
  Total_Viscosity_j = Laminar_Viscosity_j+Eddy_Viscosity_j;
  
  visc_div_vel = 0.0;
  if (!limiter)
    for (iDim = 0 ; iDim < nDim; iDim++)
      visc_div_vel += .5*(Total_Viscosity_i*PrimVar_Grad_i[iDim+1][iDim] + 
			  Total_Viscosity_j*PrimVar_Grad_j[iDim+1][iDim]);
  else
    for (iDim = 0 ; iDim < nDim; iDim++)
      visc_div_vel += .5*(Total_Viscosity_i*PrimVar_Grad_i[iDim+1][iDim]*PrimVar_Lim_i[iVar] + 
			  Total_Viscosity_j*PrimVar_Grad_j[iDim+1][iDim]*PrimVar_Lim_j[iVar]);

  for (iVar = 0; iVar < nDim; iVar++) {
    if (!limiter) {
      for (iDim = 0; iDim < nDim; iDim++) {
	Mean_Tau[iVar][iDim] = 0.5*(Total_Viscosity_i*PrimVar_Grad_i[iVar+1][iDim] + 
				    Total_Viscosity_j*PrimVar_Grad_j[iVar+1][iDim] +
				    Total_Viscosity_i*PrimVar_Grad_i[iDim+1][iVar] + 
				    Total_Viscosity_j*PrimVar_Grad_j[iDim+1][iVar]
				    );
      }
    }
    else {
      for (iDim = 0; iDim < nDim; iDim++) {
	Mean_Tau[iVar][iDim] = 0.5*(Total_Viscosity_i*PrimVar_Grad_i[iVar+1][iDim]*PrimVar_Lim_i[iVar] +
				    Total_Viscosity_j*PrimVar_Grad_j[iVar+1][iDim]*PrimVar_Lim_j[iVar] +
				    Total_Viscosity_i*PrimVar_Grad_i[iDim+1][iVar]*PrimVar_Lim_i[iDim] +
				    Total_Viscosity_j*PrimVar_Grad_j[iDim+1][iVar]*PrimVar_Lim_j[iDim]
				    );
      }
    }
  }

  for (iDim = 0 ; iDim < nDim; iDim++)
    Mean_Tau[iDim][iDim] -= TWO3*(visc_div_vel + .5*Mean_turb_ke);

  /*--- Get projected flux tensor ---*/
  
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_Tau, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
  /*--- Save residual value ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  
  if (implicit) {
    
    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
	for (jVar = 0; jVar < nVar; jVar++) {
	  val_Jacobian_i[iVar][jVar] = 0.0;
	  val_Jacobian_j[iVar][jVar] = 0.0;
	}
      }
    }
    else {
      GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
			 sqrt(dist_ij_2), UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}

//CGeneralAvgGradCorrected_Flow::CGeneralAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
//
//  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
//
//  /*--- Compressible flow, primitive variables nDim+3, (T, vx, vy, vz, P, rho) ---*/
//  PrimVar_i = new double [nDim+3];
//  PrimVar_j = new double [nDim+3];
//  Mean_PrimVar = new double [nDim+3];
//  Mean_SecVar = new double [8];
//
//  /*--- Compressible flow, primitive gradient variables nDim+1, (T, vx, vy, vz) ---*/
//  Proj_Mean_GradPrimVar_Edge = new double [nDim+1];
//  Mean_GradPrimVar = new double* [nDim+1];
//  for (iVar = 0; iVar < nDim+1; iVar++)
//    Mean_GradPrimVar[iVar] = new double [nDim];
//
//  Edge_Vector = new double [nDim];
//
//}
//
//CGeneralAvgGradCorrected_Flow::~CGeneralAvgGradCorrected_Flow(void) {
//
//  delete [] PrimVar_i;
//  delete [] PrimVar_j;
//  delete [] Mean_PrimVar;
//  delete [] Mean_SecVar;
//  delete [] Proj_Mean_GradPrimVar_Edge;
//  delete [] Edge_Vector;
//
//  for (iVar = 0; iVar < nDim+1; iVar++)
//    delete [] Mean_GradPrimVar[iVar];
//  delete [] Mean_GradPrimVar;
//
//}
//
//void CGeneralAvgGradCorrected_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
//
//  /*--- Normalized normal vector ---*/
//
//  Area = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++)
//    Area += Normal[iDim]*Normal[iDim];
//  Area = sqrt(Area);
//
//  for (iDim = 0; iDim < nDim; iDim++)
//    UnitNormal[iDim] = Normal[iDim]/Area;
//
//  /*--- Compute vector going from iPoint to jPoint ---*/
//
//  dist_ij_2 = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
//    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
//  }
//
//  /*--- Laminar and Eddy viscosity ---*/
//
//  Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
//  Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
//  Thermal_Conductivity_i = V_i[nDim+7]; Thermal_Conductivity_j = V_j[nDim+7];
//  Cp_i = V_i[nDim+8]; Cp_j = V_j[nDim+8];
//
//  for (iVar = 0; iVar < nDim+3; iVar++) {
//    PrimVar_i[iVar] = V_i[iVar];
//    PrimVar_j[iVar] = V_j[iVar];
//    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
//  }
//
//  /*--- Secondary variables ---*/
//  for (iVar = 0; iVar < 8; iVar++) {
//    Mean_SecVar[iVar] = 0.5*(S_i[iVar]+S_j[iVar]);
//  }
//
//  /*--- Mean Viscosities and turbulent kinetic energy ---*/
//
//  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
//  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
//  Mean_turb_ke              = 0.5*(turb_ke_i + turb_ke_j);
//  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);
//  Mean_Cp                   = 0.5*(Cp_i + Cp_j);
//
//  /*--- Projection of the mean gradient in the direction of the edge ---*/
//
//  for (iVar = 0; iVar < nDim+1; iVar++) {
//    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
//    for (iDim = 0; iDim < nDim; iDim++) {
//      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
//      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
//    }
//    if (dist_ij_2 != 0.0) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
//                                         (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
//      }
//    }
//  }
//
//  /*--- Get projected flux tensor ---*/
//
//  GetViscousProjFlux( Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
//		              Mean_Thermal_Conductivity, Mean_Cp );
//
//  /*--- Save residual value ---*/
//
//  for (iVar = 0; iVar < nVar; iVar++)
//    val_residual[iVar] = Proj_Flux_Tensor[iVar];
//
//  /*--- Compute the implicit part ---*/
//
//  if (implicit) {
//
//    if (dist_ij_2 == 0.0) {
//      for (iVar = 0; iVar < nVar; iVar++) {
//        for (jVar = 0; jVar < nVar; jVar++) {
//          val_Jacobian_i[iVar][jVar] = 0.0;
//          val_Jacobian_j[iVar][jVar] = 0.0;
//        }
//      }
//    }
//    else {
////		GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
////				sqrt(dist_ij_2), UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
//        GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar, Mean_SecVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_Thermal_Conductivity, Mean_Cp,
//        				sqrt(dist_ij_2), UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
//    }
//
//  }
//
//}

CGeneralAvgGradCorrected_Flow::CGeneralAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Compressible flow, primitive variables nDim+3, (vx, vy, vz, P, rho, h) ---*/
  PrimVar_i = new double [nDim+4];
  PrimVar_j = new double [nDim+4];
  Mean_PrimVar = new double [nDim+4];
  Mean_SecVar = new double [2];
  
  /*--- Compressible flow, primitive gradient variables nDim+1, (T, vx, vy, vz) ---*/
  Proj_Mean_GradPrimVar_Edge = new double [nDim+1];
  Mean_GradPrimVar = new double* [nDim+1];
  for (iVar = 0; iVar < nDim+1; iVar++)
    Mean_GradPrimVar[iVar] = new double [nDim];
  
  Edge_Vector = new double [nDim];
  
}

CGeneralAvgGradCorrected_Flow::~CGeneralAvgGradCorrected_Flow(void) {
  
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  delete [] Mean_SecVar;
  delete [] Proj_Mean_GradPrimVar_Edge;
  delete [] Edge_Vector;
  
  for (iVar = 0; iVar < nDim+1; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
}

void CGeneralAvgGradCorrected_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  Thermal_Conductivity_i = V_i[nDim+7]; Thermal_Conductivity_j = V_j[nDim+7];
  Cp_i = V_i[nDim+8]; Cp_j = V_j[nDim+8];
  
  for (iVar = 0; iVar < nDim+4; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  
  /*--- Secondary variables ---*/
  for (iVar = 0; iVar < 2; iVar++) {
    Mean_SecVar[iVar] = 0.5*(S_i[iVar+2]+S_j[iVar+2]);
  }

  /*--- Mean Viscosities and turbulent kinetic energy ---*/
  
  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke              = 0.5*(turb_ke_i + turb_ke_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);
  Mean_Cp                   = 0.5*(Cp_i + Cp_j);
  
  /*--- Projection of the mean gradient in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nDim+1; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    }
    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                         (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
      }
    }
  }
  
  /*--- Get projected flux tensor ---*/
  
  GetViscousProjFlux( Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
		              Mean_Thermal_Conductivity, Mean_Cp );
  
  /*--- Save residual value ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  
  if (implicit) {
    
    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
//		GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
//				sqrt(dist_ij_2), UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
        GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar, Mean_SecVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_Thermal_Conductivity, Mean_Cp,
        				sqrt(dist_ij_2), UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}

CAvgGradCorrectedArtComp_Flow::CAvgGradCorrectedArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  PrimVar_i = new double [nVar];
  PrimVar_j = new double [nVar];
  Proj_Mean_GradPrimVar_Edge = new double [nVar];
  Edge_Vector = new double [nDim];
  
  Mean_GradPrimVar = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPrimVar[iVar] = new double [nDim];
  
}

CAvgGradCorrectedArtComp_Flow::~CAvgGradCorrectedArtComp_Flow(void) {
  
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Proj_Mean_GradPrimVar_Edge;
  delete [] Edge_Vector;
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
}

void CAvgGradCorrectedArtComp_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Conversion to Primitive Variables ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
  Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  
  /*--- Mean Viscosities ---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Projection of the mean gradient in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    }
    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                         (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
      }
    }
  }
  
  /*--- Get projected flux tensor ---*/
  
  GetViscousArtCompProjFlux(Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Implicit part ---*/
  
  if (implicit) {
    
    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, sqrt(dist_ij_2), UnitNormal,
                                Area, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}

CSourceGravity::CSourceGravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  freesurface = (config->GetKind_Regime() == FREESURFACE);
  
}

CSourceGravity::~CSourceGravity(void) { }

void CSourceGravity::ComputeResidual(double *val_residual, CConfig *config) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = 0.0;
  
  if (compressible) {
    
    /*--- Evaluate the source term  ---*/
    val_residual[nDim] = Volume * U_i[0] * STANDART_GRAVITY;
    
  }
  if (incompressible || freesurface) {
    
    /*--- Compute the Froude number  ---*/
    Froude = config->GetFroude();
    
    /*--- Evaluate the source term  ---*/
    val_residual[nDim] = Volume * DensityInc_i / (Froude * Froude);
    
  }
  
}

CSourceRotatingFrame_Flow::CSourceRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
}

CSourceRotatingFrame_Flow::~CSourceRotatingFrame_Flow(void) { }

void CSourceRotatingFrame_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
  
  unsigned short iDim, iVar, jVar;
  double Omega[3] = {0,0,0}, Momentum[3] = {0,0,0};
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Retrieve the angular velocity vector from config. ---*/
  
  Omega[0]  = config->GetRotation_Rate_X(ZONE_0)/config->GetOmega_Ref();
  Omega[1]  = config->GetRotation_Rate_Y(ZONE_0)/config->GetOmega_Ref();
  Omega[2]  = config->GetRotation_Rate_Z(ZONE_0)/config->GetOmega_Ref();
  
  /*--- Get the momentum vector at the current node. ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Momentum[iDim] = U_i[iDim+1];
  
  /*--- Calculate rotating frame source term as ( Omega X Rho-U ) ---*/
  
  if (nDim == 2) {
    val_residual[0] = 0.0;
    val_residual[1] = (Omega[1]*Momentum[2] - Omega[2]*Momentum[1])*Volume;
    val_residual[2] = (Omega[2]*Momentum[0] - Omega[0]*Momentum[2])*Volume;
    val_residual[3] = 0.0;
  } else {
    val_residual[0] = 0.0;
    val_residual[1] = (Omega[1]*Momentum[2] - Omega[2]*Momentum[1])*Volume;
    val_residual[2] = (Omega[2]*Momentum[0] - Omega[0]*Momentum[2])*Volume;
    val_residual[3] = (Omega[0]*Momentum[1] - Omega[1]*Momentum[0])*Volume;
    val_residual[4] = 0.0;
  }
  
  /*--- Calculate the source term Jacobian ---*/
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
    if (nDim == 2) {
      val_Jacobian_i[1][2] = -Omega[2]*Volume;
      val_Jacobian_i[2][1] =  Omega[2]*Volume;
    } else {
      val_Jacobian_i[1][2] = -Omega[2]*Volume;
      val_Jacobian_i[1][3] =  Omega[1]*Volume;
      val_Jacobian_i[2][1] =  Omega[2]*Volume;
      val_Jacobian_i[2][3] = -Omega[0]*Volume;
      val_Jacobian_i[3][1] = -Omega[1]*Volume;
      val_Jacobian_i[3][2] =  Omega[0]*Volume;
    }
  }
  
}

CSourceAxisymmetric_Flow::CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
}

CSourceAxisymmetric_Flow::~CSourceAxisymmetric_Flow(void) { }

void CSourceAxisymmetric_Flow::ComputeResidual(double *val_residual, double **Jacobian_i, CConfig *config) {
  
  double yinv, Pressure_i, Enthalpy_i, Velocity_i, sq_vel;
  unsigned short iDim;
  
  bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  
  if (Coord_i[1] > 0.0) yinv = 1.0/Coord_i[1];
  else yinv = 0.0;
  
  if (compressible) {
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i = U_i[iDim+1] / U_i[0];
      sq_vel += Velocity_i *Velocity_i;
    }
    
    Pressure_i = (Gamma-1.0)*U_i[0]*(U_i[nDim+1]/U_i[0]-0.5*sq_vel);
    Enthalpy_i = (U_i[nDim+1] + Pressure_i) / U_i[0];
    
    val_residual[0] = yinv*Volume*U_i[2];
    val_residual[1] = yinv*Volume*U_i[1]*U_i[2]/U_i[0];
    val_residual[2] = yinv*Volume*(U_i[2]*U_i[2]/U_i[0]);
    val_residual[3] = yinv*Volume*Enthalpy_i*U_i[2];
  }
  if (incompressible || freesurface) {
    val_residual[0] = yinv*Volume*U_i[2]*BetaInc2_i;
    val_residual[1] = yinv*Volume*U_i[1]*U_i[2]/DensityInc_i;
    val_residual[2] = yinv*Volume*U_i[2]*U_i[2]/DensityInc_i;
  }
  
  
  if (implicit) {
    Jacobian_i[0][0] = 0;
    Jacobian_i[0][1] = 0;
    Jacobian_i[0][2] = 1.;
    Jacobian_i[0][3] = 0;
    
    Jacobian_i[1][0] = -U_i[1]*U_i[2]/(U_i[0]*U_i[0]);
    Jacobian_i[1][1] = U_i[2]/U_i[0];
    Jacobian_i[1][2] = U_i[1]/U_i[0];
    Jacobian_i[1][3] = 0;
    
    Jacobian_i[2][0] = -U_i[2]*U_i[2]/(U_i[0]*U_i[0]);
    Jacobian_i[2][1] = 0;
    Jacobian_i[2][2] = 2*U_i[2]/U_i[0];
    Jacobian_i[2][3] = 0;
    
    Jacobian_i[3][0] = -Gamma*U_i[2]*U_i[3]/(U_i[0]*U_i[0]) + (Gamma-1)*U_i[2]*(U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]*U_i[0]);
    Jacobian_i[3][1] = -(Gamma-1)*U_i[2]*U_i[1]/(U_i[0]*U_i[0]);
    Jacobian_i[3][2] = Gamma*U_i[3]/U_i[0] - 1/2*(Gamma-1)*( (U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]) + 2*U_i[2]*U_i[2]/(U_i[0]*U_i[0]) );
    Jacobian_i[3][3] = Gamma*U_i[2]/U_i[0];
    
    for (int iVar=0; iVar<4; iVar++)
      for (int jVar=0; jVar<4; jVar++)
        Jacobian_i[iVar][jVar] *= yinv*Volume;
  }
  
}

CSourceWindGust::CSourceWindGust(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
}

CSourceWindGust::~CSourceWindGust(void) { }

void CSourceWindGust::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) {
  
  double u_gust, v_gust, du_gust_dx, du_gust_dy, du_gust_dt, dv_gust_dx, dv_gust_dy, dv_gust_dt, smx, smy, se, rho, u, v, p;
  unsigned short GustDir = config->GetGust_Dir(); //Gust direction
  
  u_gust = WindGust_i[0];
  v_gust = WindGust_i[1];
  
  if (GustDir == X_DIR) {
    du_gust_dx = WindGustDer_i[0];
    du_gust_dy = WindGustDer_i[1];
    du_gust_dt = WindGustDer_i[2];
    dv_gust_dx = 0.0;
    dv_gust_dy = 0.0;
    dv_gust_dt = 0.0;
  } else {
    du_gust_dx = 0.0;
    du_gust_dy = 0.0;
    du_gust_dt = 0.0;
    dv_gust_dx = WindGustDer_i[0];
    dv_gust_dy = WindGustDer_i[1];
    dv_gust_dt = WindGustDer_i[2];
    
  }
  
  /*--- Primitive variables at point i ---*/
  u = V_i[1];
  v = V_i[2];
  p = V_i[nDim+1];
  rho = V_i[nDim+2];
  
  /*--- Source terms ---*/
  smx = rho*(du_gust_dt + (u+u_gust)*du_gust_dx + (v+v_gust)*du_gust_dy);
  smy = rho*(dv_gust_dt + (u+u_gust)*dv_gust_dx + (v+v_gust)*dv_gust_dy);
  se = u*smx + v*smy + p*(du_gust_dx + dv_gust_dy);
  
  if (nDim == 2) {
    val_residual[0] = 0.0;
    val_residual[1] = smx*Volume;
    val_residual[2] = smy*Volume;
    val_residual[3] = se*Volume;
  } else {
    cout << "ERROR: You should only be in the gust source term in two dimensions" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
	MPI_Abort(MPI_COMM_WORLD,1);
	MPI_Finalize();
#endif
    
  }
  
  /*--- For now the source term Jacobian is just set to zero ---*/
  
  unsigned short iVar, jVar;
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Calculate the source term Jacobian ---*/
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
  }
  
}


