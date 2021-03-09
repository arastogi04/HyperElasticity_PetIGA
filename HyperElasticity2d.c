#include "petiga.h"
#include "petscmat.h"
#include "petscblaslapack.h"
#include "slepceps.h"
/*
  This code implements a HyperElastic material model in the context of
  large deformation elasticity. Implementation credit goes to students
  of the 2012 summer course `Nonlinear Finite Element Analysis' given
  in Universidad de los Andes, Bogota, Colombia:

  Lina Mar¿a Bernal Martinez
  Gabriel Andres Espinosa Barrios
  Federico Fuentes Caycedo
  Juan Camilo Mahecha Zambrano

 */

typedef struct {
  PetscReal lambda,mu,a,b,c,d,c1,c2,kappa;
  void (*model) (PetscScalar u[2], PetscScalar grad_u[2][2], PetscScalar (*F)[2], PetscScalar (*S)[2], PetscScalar (*D)[3], void *ctx);
} AppCtx;

static void NeoHookeanModel(PetscScalar u[2], PetscScalar grad_u[2][2], PetscScalar (*F)[2], PetscScalar (*S)[2], PetscScalar (*D)[3], void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

  PetscReal lambda = user->lambda;
  PetscReal mu = user->mu;
  // F = I + u_{i,j}
  F[0][0] = 1+grad_u[0][0]; F[0][1] =   grad_u[0][1];  
  F[1][0] =   grad_u[1][0]; F[1][1] = 1+grad_u[1][1];  
  

  // Finv
  PetscScalar Finv[2][2],J,Jinv;
  J  = F[0][0]*F[1][1] - F[1][0]*F[0][1];
  Jinv = 1./J;
  Finv[0][0] =  (F[1][1])*Jinv;
  Finv[1][0] = -(F[1][0])*Jinv;
  Finv[0][1] = -(F[0][1])*Jinv;
  Finv[1][1] =  (F[0][0])*Jinv;

  // C^-1 = (F^T F)^-1 = F^-1 F^-T
  PetscScalar Cinv[3][3];
  Cinv[0][0] = Finv[0][0]*Finv[0][0] + Finv[0][1]*Finv[0][1];
  Cinv[0][1] = Finv[0][0]*Finv[1][0] + Finv[0][1]*Finv[1][1];
  Cinv[1][0] = Finv[1][0]*Finv[0][0] + Finv[1][1]*Finv[0][1];
  Cinv[1][1] = Finv[1][0]*Finv[1][0] + Finv[1][1]*Finv[1][1];
  
// 2nd Piola-Kirchoff stress tensor
  PetscScalar temp=(0.5*lambda)*(J*J-1.0);
  S[0][0] = temp*Cinv[0][0] + mu*(1.0-Cinv[0][0]);
  S[0][1] = temp*Cinv[0][1] + mu*(-Cinv[0][1]);
  //S[0][2] = temp*Cinv[0][2] + mu*(-Cinv[0][2]);
  S[1][0] = temp*Cinv[1][0] + mu*(-Cinv[1][0]);
  S[1][1] = temp*Cinv[1][1] + mu*(1.0-Cinv[1][1]);
  //S[1][2] = temp*Cinv[1][2] + mu*(-Cinv[1][2]);
  //S[2][0] = temp*Cinv[2][0] + mu*(-Cinv[2][0]);
 // S[2][1] = temp*Cinv[2][1] + mu*(-Cinv[2][1]);
 // S[2][2] = temp*Cinv[2][2] + mu*(1.0-Cinv[2][2]);

  // C_abcd=lambda*J^2*Cinv_ab*Cinv_cd+[2*miu-lambda(J^2-1)]*0.5(Cinv_ac*Cinv_bd+Cinv_ad*Cinv_bc)
  PetscScalar temp1=lambda*J*J;
  PetscScalar temp2=2*mu-lambda*(J*J-1);
  D[0][0]=temp1*Cinv[0][0]*Cinv[0][0]+temp2*0.5*(Cinv[0][0]*Cinv[0][0]+Cinv[0][0]*Cinv[0][0]);
  D[0][1]=temp1*Cinv[0][0]*Cinv[1][1]+temp2*0.5*(Cinv[0][1]*Cinv[0][1]+Cinv[0][1]*Cinv[0][1]);
  //D[0][2]=temp1*Cinv[0][0]*Cinv[2][2]+temp2*0.5*(Cinv[0][2]*Cinv[0][2]+Cinv[0][2]*Cinv[0][2]);
  D[0][2]=temp1*Cinv[0][0]*Cinv[0][1]+temp2*0.5*(Cinv[0][0]*Cinv[0][1]+Cinv[0][1]*Cinv[0][0]);
  D[1][1]=temp1*Cinv[1][1]*Cinv[1][1]+temp2*0.5*(Cinv[1][1]*Cinv[1][1]+Cinv[1][1]*Cinv[1][1]);
  D[1][2]=temp1*Cinv[1][1]*Cinv[0][1]+temp2*0.5*(Cinv[1][0]*Cinv[1][1]+Cinv[1][1]*Cinv[1][0]);
  D[2][2]=temp1*Cinv[0][1]*Cinv[0][1]+temp2*0.5*(Cinv[0][0]*Cinv[1][1]+Cinv[0][1]*Cinv[1][0]);
  D[1][0]=D[0][1];
  D[2][0]=D[0][2];
  D[2][1]=D[1][2];
}


static void DeltaE(PetscScalar Nx, PetscScalar Ny, PetscScalar (*F)[2], PetscScalar (*B)[2])
{
  // Given F and basis values, returns B
  B[0][0] = F[0][0]*Nx;
  B[0][1] = F[1][0]*Nx;
  B[1][0] = F[0][1]*Ny;
  B[1][1] = F[1][1]*Ny;
  B[2][0] = F[0][0]*Ny+F[0][1]*Nx;
  B[2][1] = F[1][0]*Ny+F[1][1]*Nx;

}


static PetscErrorCode Residual(IGAPoint pnt,const PetscScalar *U,PetscScalar *Re,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

  // call user model
  PetscScalar F[2][2],S[2][2],D[3][3],B[3][2];

  // Interpolate the solution and gradient given U
  PetscScalar u[2];
  PetscScalar grad_u[2][2];
  IGAPointFormValue(pnt,U,&u[0]);
  IGAPointFormGrad (pnt,U,&grad_u[0][0]);


  user->model(u,grad_u,F,S,D,ctx);

  // Get basis function gradients
  PetscReal (*N1)[2] = (PetscReal (*)[2]) pnt->shape[1];

  PetscScalar (*R)[2] = (PetscScalar (*)[2])Re;
  PetscInt a,nen=pnt->nen;
  for (a=0; a<nen; a++) {
    PetscReal Na_x  = N1[a][0];
    PetscReal Na_y  = N1[a][1];
    DeltaE(Na_x,Na_y,F,B);
    R[a][0] = B[0][0]*S[0][0]+B[1][0]*S[1][1];
    R[a][1] = B[0][1]*S[0][0]+B[1][1]*S[1][1];
  }
  return 0;
}

static PetscErrorCode Jacobian(IGAPoint pnt,const PetscScalar *U,PetscScalar *Je,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

  // Call user model
  PetscScalar F[2][2],S[2][2],D[3][3];
  
  // Interpolate the solution and gradient given U
  PetscScalar u[2];
  PetscScalar grad_u[2][2];
  IGAPointFormValue(pnt,U,&u[0]);
  IGAPointFormGrad (pnt,U,&grad_u[0][0]);


  user->model(u,grad_u,F,S,D,ctx);

  // Get basis function gradients
  PetscReal (*N1)[2] = (PetscReal (*)[2]) pnt->shape[1];

  // Put together the jacobian
  PetscInt a,b,nen=pnt->nen;
  PetscScalar (*K)[2][nen][2] = (PetscScalar (*)[2][nen][2])Je;
  PetscScalar G;
  PetscScalar Chi[3][3];
  PetscScalar Ba[3][2];
  PetscScalar Bb[3][2];

  for (a=0; a<nen; a++) {

    PetscReal Na_x  = N1[a][0];
    PetscReal Na_y  = N1[a][1];
    DeltaE(Na_x,Na_y,F,Ba);

    for (b=0; b<nen; b++) {

      PetscReal Nb_x  = N1[b][0];
      PetscReal Nb_y  = N1[b][1];
      DeltaE(Nb_x,Nb_y,F,Bb);

      Chi[0][0]=Bb[0][0]*(Ba[0][0]*D[0][0] + Ba[1][0]*D[1][0] + Ba[2][0]*D[2][0]) +
        Bb[1][0]*(Ba[0][0]*D[0][1] + Ba[1][0]*D[1][1] + Ba[2][0]*D[2][1]) +
        Bb[2][0]*(Ba[0][0]*D[0][2] + Ba[1][0]*D[1][2] + Ba[2][0]*D[2][2]);

      Chi[0][1]=Bb[0][1]*(Ba[0][0]*D[0][0] + Ba[1][0]*D[1][0] + Ba[2][0]*D[2][0]) +
        Bb[1][1]*(Ba[0][0]*D[0][1] + Ba[1][0]*D[1][1] + Ba[2][0]*D[2][1]) +
        Bb[2][1]*(Ba[0][0]*D[0][2] + Ba[1][0]*D[1][2] + Ba[2][0]*D[2][2]);


      Chi[1][0]=Bb[0][0]*(Ba[0][1]*D[0][0] + Ba[1][1]*D[1][0] + Ba[2][1]*D[2][0]) +
        Bb[1][0]*(Ba[0][1]*D[0][1] + Ba[1][1]*D[1][1] + Ba[2][1]*D[2][1] ) +
        Bb[2][0]*(Ba[0][1]*D[0][2] + Ba[1][1]*D[1][2] + Ba[2][1]*D[2][2]);

      Chi[1][1]=Bb[0][1]*(Ba[0][1]*D[0][0] + Ba[1][1]*D[1][0] + Ba[2][1]*D[2][0]) +
        Bb[1][1]*(Ba[0][1]*D[0][1] + Ba[1][1]*D[1][1] + Ba[2][1]*D[2][1] ) +
        Bb[2][1]*(Ba[0][1]*D[0][2] + Ba[1][1]*D[1][2] + Ba[2][1]*D[2][2] ) ;

             G=Na_x*(S[0][0]*Nb_x + S[0][1]*Nb_y) +
        Na_y*(S[1][0]*Nb_x + S[1][1]*Nb_y);

      K[a][0][b][0] = G+Chi[0][0];
      K[a][1][b][0] =   Chi[1][0];

      K[a][0][b][1] =   Chi[0][1];
      K[a][1][b][1] = G+Chi[1][1];

    }
  }

  return 0;
}

int main(int argc, char *argv[])
{
  // Initialization of PETSc
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  // Application specific data (defaults to Aluminum)
  AppCtx user;
  PetscScalar E  = 10;
  PetscScalar nu = 0.3;
  PetscInt nsteps = 1;

 
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","HyperElasticity Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nsteps","Number of load steps to take",__FILE__,nsteps,&nsteps,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  user.model  = NeoHookeanModel;
  user.lambda = E*nu/(1+nu)/(1-2*nu);
  user.mu     = 0.5*E/(1+nu);


  // Initialize the discretization
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDof(iga,2);CHKERRQ(ierr);
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);

  if(iga->geometry == 0 && nsteps > 50){
    SETERRQ(PETSC_COMM_WORLD,
            PETSC_ERR_ARG_OUTOFRANGE,
            "You must specify a geometry to use an updated Lagrangian approach");
  }

  // Set boundary conditions
  ierr = IGASetBoundaryValue(iga,0,0,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,1,0,-0.1/((PetscReal)nsteps));CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,1,0,1,0.0);CHKERRQ(ierr);
  //ierr = IGASetBoundaryValue(iga,2,0,2,0.0);CHKERRQ(ierr);
  //ierr = IGASetBoundaryValue(iga,2,1,2,0.0);CHKERRQ(ierr);
  

  // Setup the nonlinear solver
  SNES snes;
  KSP ksp;
  PC pc;
  Vec U,Utotal;
  ierr = IGASetFormFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormJacobian(iga,Jacobian,&user);CHKERRQ(ierr);
  ierr = IGACreateSNES(iga,&snes);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&Utotal);CHKERRQ(ierr);
  ierr = VecZeroEntries(Utotal);CHKERRQ(ierr);
  ierr = IGAWrite(iga,"geometry0.dat");CHKERRQ(ierr);
  ierr = IGAWriteVec(iga,Utotal,"disp0.dat");CHKERRQ(ierr);
  ierr = SNESGetKSP(snes, &ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
   
    

  // Load stepping
  PetscInt step;
 
  for(step=0;step<nsteps;step++){

    PetscPrintf(PETSC_COMM_WORLD,"%d Load Step\n",step);
    Mat J;

    // Solve step
    ierr = VecZeroEntries(U);CHKERRQ(ierr);
    ierr = SNESSolve(snes,NULL,U);CHKERRQ(ierr);

    //Get the Consistent Tangent
     SNESGetJacobian(snes, &J, NULL, NULL, NULL);
     
    // Store total displacement
     ierr = VecAXPY(Utotal,1.0,U);CHKERRQ(ierr);
     PetscViewer viewer;
     PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Amat.m", &viewer);
     PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
     MatView(J,viewer);
     PetscViewerPopFormat(viewer);
     PetscViewerDestroy(&viewer);
    
    // Update the geometry
    if(iga->geometry){
      Vec localU;
      const PetscScalar *arrayU;
      ierr = IGAGetLocalVecArray(iga,U,&localU,&arrayU);CHKERRQ(ierr);
      PetscInt i,N;
      ierr = VecGetSize(localU,&N);
      for(i=0;i<N;i++) iga->geometryX[i] += arrayU[i];
      ierr = IGARestoreLocalVecArray(iga,U,&localU,&arrayU);CHKERRQ(ierr);
    }

    // Dump solution vector
    char filename[256];
    sprintf(filename,"disp%d.dat",step+1);
    ierr = IGAWriteVec(iga,Utotal,filename);CHKERRQ(ierr);

    //Write the geometry
    char filenamegeo[256];
    sprintf(filenamegeo,"geometry%d.dat",step+1);
    ierr = IGAWrite(iga,filenamegeo);CHKERRQ(ierr);
      
}

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
