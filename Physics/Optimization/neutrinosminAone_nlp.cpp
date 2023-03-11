// neutrinos.cpp
// Nonlinear Ipopt program

// Author: Bryan A. Toth
// btoth@physics.ucsd.edu

#include "neutrinosminAone_nlp.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cstring>
#ifndef HAVE_CSTDIO
#define HAVE_CSTDIO
# include <cstdio>
# include <iostream>
# include <fstream>
# include <string>
# include <stdlib.h>
#else
# ifndef HAVE_STDIO_H
#define HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

using namespace Ipopt;

using namespace std;

// constructor
NEUTRINOS_NLP::NEUTRINOS_NLP(int id)
{
  nU=0;
  nP=0;
  nY=2;
  nM=3;
  nI=0;

  K11val = new double[nU];
  K11val2 = new double[nU];
  K11valp1 = new double[nU];
  dK11val = new double[nU];
  dK11val2 = new double[nU];
  dK11valp1 = new double[nU];
  Xdval = new double[nM];
  Xdval2 = new double[nM];
  Xdvalp1 = new double[nM];
  Xval = new double[nY];
  Xval2 = new double[nY];
  Xvalp1 = new double[nY];
  Pval = new double[nP];
  Ival = new double[nI];
  Ival2 = new double[nI];
  Ivalp1 = new double[nI];
  Rf0 = new double[nY];

  string buffer;
  specs = new string[9+nP+nY+nI+nU+nM+1+1];

  int count;
  count = 0;
  
  ifstream fin ("specs.txt");
  if (fin.is_open()){
    while (! fin.eof()){
      getline (fin,buffer);
      if (buffer[0] !='#'){
        specs[count] = buffer;
        count++;
      }
    }
    fin.close();
  }
  else cout << "Unable to open file";

  Time = atoi(specs[0].c_str());
  skip = atoi(specs[1].c_str());
  hstep = atof(specs[2].c_str());

  string filename;
  int ret;
  double tempdata;
  Ntotal = (2*Time+1)*nY +nP;
  solution = new double[Ntotal];
  int toggle = 0;
  file_fmt = atoi(specs[3].c_str());

  // Read in observation files

    if (file_fmt == 1) toggle = 1;
    if (file_fmt == 2) toggle = 2;
    if (file_fmt == 3) toggle = 3;

    if (file_fmt < 2) {
      // Read in files one by one if file_fmt = 0,1

      taskid = id;

      dd = new double[2*Time+1];
      dddummy = new double[skip];
      FILE *pFile0;
      filename = specs[4+toggle];
      pFile0 = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++)
      {
        ret = fscanf (pFile0, "%lf", &dddummy[jt]);
        if (ret == EOF) break;
      }
      for(Index jt=0;jt<2*Time+1;jt++)
      {
        ret = fscanf (pFile0, "%lf", &dd[jt]);
        if (ret == EOF) break;
      }
      fclose (pFile0);

      BC = new double[2*Time+1];
      BCdummy = new double[skip];
      FILE *pFile1;
      filename = specs[5+toggle];
      pFile1 = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++)
      {
        ret = fscanf (pFile1, "%lf", &BCdummy[jt]);
        if (ret == EOF) break;
      }
      for(Index jt=0;jt<2*Time+1;jt++)
      {
        ret = fscanf (pFile1, "%lf", &BC[jt]);
        if (ret == EOF) break;
      }
      fclose (pFile1);

      tt = new double[2*Time+1];
      ttdummy = new double[skip];
      FILE *pFile2;
      filename = specs[6+toggle];
      pFile2 = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++)
      {
        ret = fscanf (pFile2, "%lf", &ttdummy[jt]);
        if (ret == EOF) break;
      }
      for(Index jt=0;jt<2*Time+1;jt++)
      {
        ret = fscanf (pFile2, "%lf", &tt[jt]);
        if (ret == EOF) break;
      }
      fclose (pFile2);

    }
    else if (file_fmt >= 2){
      // Read in files from single file if file_fmt = 2,3

      dd = new double[2*Time+1];
      dddummy = new double[skip];
      BC = new double[2*Time+1];
      BCdummy = new double[skip];
      tt = new double[2*Time+1];
      ttdummy = new double[skip];

      FILE *pFile0;
      filename = specs[4+toggle];
      char init_idx[21]; 
      modid = atoi(specs[4].c_str());
      if (modid != 0){
        pathid = id / modid;
        taskid = id % modid;
        sprintf(init_idx, "%d.%s", pathid, specs[5].c_str());
      }
      else{
        pathid = 0;
        taskid = id;
        sprintf(init_idx, ".%s", specs[5].c_str());
      }
      filename += init_idx;
      pFile0 = fopen(filename.c_str(),"r");

      for(Index rows=0;rows<skip;rows++){
        for (Index cols=0;cols<nY;cols++){
          if (cols == 0){
            ret = fscanf (pFile0, "%lf", &dddummy[rows]);
          }
          else if (cols == 1){
            ret = fscanf (pFile0, "%lf", &BCdummy[rows]);
          }
          else if (cols == 2){
            ret = fscanf (pFile0, "%lf", &ttdummy[rows]);
          }
          else {
            ret = fscanf (pFile0, "%lf", &tempdata);
          }          
           if (ret == EOF) break;
        }
      }
      for(Index rows=0;rows<2*Time+1;rows++){
        for (Index cols=0;cols<nY;cols++){
          if (cols == 0){
            ret = fscanf (pFile0, "%lf", &dd[rows]);
          }
          else if (cols == 1){
            ret = fscanf (pFile0, "%lf", &BC[rows]);
          }
          else if (cols == 2){
            ret = fscanf (pFile0, "%lf", &tt[rows]);
          }
          else {
            ret = fscanf (pFile0, "%lf", &tempdata);
          }          
           if (ret == EOF) break;
        }
      }
      fclose (pFile0);
    }

  //Read in boundary conditions

  if (file_fmt == 0) toggle = nM + nI;
  if (file_fmt == 1) toggle = nM + nI + 1;
  if (file_fmt == 2) toggle = (nM > 0) + (nI > 0) + 2;
  if (file_fmt == 3) toggle = (nM > 0) + (nI > 0) + 3;

  int rows = nY+nU+nP+1;
  bounds = new double*[rows];
  for (Index i=0;i<rows;i++) bounds[i] = new double[4];
  int counter;
  for(Index k=0;k<rows;k++){
    counter=0;
    char* tmp = new char[specs[4+toggle+k].size()+1];
    strcpy( tmp, specs[4+toggle+k].c_str() );
    char *ptr = strtok(tmp,",");
    bounds[k][3] = 0.0;
    while(ptr != 0){
      if(counter<3){
        bounds[k][counter] = atof(ptr);
      }
      if(counter==3) {
        bounds[k][counter] = atof(ptr);
      }
      ptr = strtok(0,",");
      counter++;
    }
  }
  for (Index i=0;i<nY;i++){
    Rf0[i]=bounds[i][2];
    bounds[i][3]=bounds[i][2];
  }

  beta=0;
  alpha = bounds[nY+nU+nP][0];
  delta_beta=(int) bounds[nY+nU+nP][1];
  max_beta=(int) bounds[nY+nU+nP][2];

  //Read in output_fmt (zero is default)
  output_fmt = atoi(specs[4+toggle+rows].c_str());
}

// destructor
NEUTRINOS_NLP::~NEUTRINOS_NLP()
{
  delete [] K11val;
  delete [] K11val2;
  delete [] K11valp1;
  delete [] dK11val;
  delete [] dK11val2;
  delete [] dK11valp1;
  delete [] Xdval;
  delete [] Xdval2;
  delete [] Xdvalp1;
  delete [] Xval;
  delete [] Xval2;
  delete [] Xvalp1;
  delete [] Pval;
  delete [] Ival;
  delete [] Ival2;
  delete [] Ivalp1;
  delete [] specs;
  delete [] dd;
  delete [] dddummy;
  delete [] BC;
  delete [] BCdummy;
  delete [] tt;
  delete [] ttdummy;
  int rows = nY+nU+nP;
  for (Index i=0;i<rows;i++) delete [] bounds[i];
  delete [] bounds;
}

bool NEUTRINOS_NLP::changeRf(){
  if((beta+delta_beta)>(max_beta-1))
    return false;
  else
    beta = beta + delta_beta;
  for (Index i=0;i<nY;i++) {
    bounds[i][3]=pow(alpha,beta)*Rf0[i];
  }
  printf("\n\n\n\n<------       Beta=%d       ------->\n\n\n\n",beta);
  return true;
}

// returns the size of the problem
bool NEUTRINOS_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag, IndexStyleEnum& index_style)

{
  // Number of variables
  n = 4*Time+2;

  // Number of equality constraints
  m = 0;

  // Number of Jacobian nonzero entries
  nnz_jac_g = 0;

  // Number of Hessian nonzero entries
  nnz_h_lag = 3*(Time+1)+14*Time+0;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}


// returns the variable bounds
bool NEUTRINOS_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // Here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  assert(n == 4*Time+2);
  assert(m == 0);

  for(Index jt=0;jt<Time+1;jt++) {
    for(Index var=0;var<nY;var++) {
      // Bounds for x
      x_l[(Time+1)*var+jt]=bounds[var][0];
      x_u[(Time+1)*var+jt]=bounds[var][1];
      // Bounds for midpoints
      if(jt<Time){
        x_l[(Time+1)*(nY+nU)+Time*var+jt]=bounds[var][0];
        x_u[(Time+1)*(nY+nU)+Time*var+jt]=bounds[var][1];
      }
    }
    for(Index cup=0;cup<nU;cup++) {
      // Bounds for k
      x_l[(Time+1)*(nY+cup)+jt]=bounds[nY+cup][0];
      x_u[(Time+1)*(nY+cup)+jt]=bounds[nY+cup][1];
      // Bounds for midpoints
      if(jt<Time) {
        x_l[(Time+1)*(nY+nU)+Time*(nY+cup)+jt]=bounds[nY+cup][0];
        x_u[(Time+1)*(nY+nU)+Time*(nY+cup)+jt]=bounds[nY+cup][1];
      }
    }
  } 

  for(Index par=0;par<nP;par++) {
    // Bounds for parameters
    x_l[2*Time*(nY+nU)+nY+nU+par]=bounds[nY+nU+par][0];
    x_u[2*Time*(nY+nU)+nY+nU+par]=bounds[nY+nU+par][1];
  }

  return true;
}
// returns the initial point for the problem
bool NEUTRINOS_NLP::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  for (Index i=0; i<n; i++) {
    x[i] = 0.0;
  }

  if(beta==0){
    double **skipinit = new double* [skip];
    double **init = new double* [2*Time + 1];
    double *param_init = new double [nP];

    for(Index i=0;i<skip;i++) skipinit[i] = new double[nY];
    for(Index i=0;i<2*Time+1;i++) init[i] = new double[nY];

    string filename;
    int ret;

    // Read in initial data files for file_fmt = 1,3

    if (file_fmt == 1){
      filename = specs[4];
      FILE *initFILE;
      initFILE = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++){
        for (Index jy=0;jy<nY;jy++){
          ret = fscanf (initFILE,"%lf",&skipinit[jt][jy]);
          if (ret == EOF) break;
        }
      }
      for(Index jt=0;jt<(2*Time+1);jt++){
        for (Index jy=0;jy<nY;jy++){
          ret = fscanf (initFILE,"%lf",&init[jt][jy]);
          if (ret == EOF) break;
        }
      }
      for(Index i=0;i<nP;i++){
        ret = fscanf (initFILE,"%lf",&param_init[i]);
        if (ret == EOF) break;
      }
      fclose (initFILE);
    }
    else if (file_fmt == 3){
      filename = specs[6];
      FILE *initFILE;
      char init_idx[21];
      if (modid != 0){
        sprintf(init_idx, "%d.%s", taskid, specs[5].c_str());
      }
      else{
        sprintf(init_idx, ".%s", specs[5].c_str());
      }
      filename += init_idx;
      initFILE = fopen(filename.c_str(),"r");

      for(Index jt=0;jt<skip;jt++){
        for (Index jy=0;jy<nY;jy++){
          ret = fscanf (initFILE,"%lf",&skipinit[jt][jy]);
          if (ret == EOF) break;
        }
      }
      for(Index jt=0;jt<(2*Time+1);jt++){
        for (Index jy=0;jy<nY;jy++){
          ret = fscanf (initFILE,"%lf",&init[jt][jy]);
          if (ret == EOF) break;
        }
      }
      for(Index i=0;i<nP;i++){
        ret = fscanf (initFILE,"%lf",&param_init[i]);
        if (ret == EOF) break;
      }
      fclose (initFILE);
    }

    // Save init_data to variables for file_fmt = 0,1,2,3

    if ((file_fmt == 1) || (file_fmt == 3)){
      for (Index jt=0;jt<Time+1;jt++){
        for (Index var=0;var<nY;var++){

        // Initial conditions for x and midpoints
          x[(Time+1)*var+jt] = init[2*jt][var];
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*var+jt] = init[2*jt+1][var];
        }

        // Initial conditions for controls and midpoints
        for (Index cup=0;cup<nU;cup++){
          x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];
        }
      }

      // Initial conditions for parameters
      for (Index par=0;par<nP;par++){
        x[(2*Time+1)*(nY+nU)+par] = param_init[par];
      }
    }
    else if ((file_fmt == 0) || (file_fmt == 2)){
      srand(taskid);
      for (Index jt=0;jt<Time+1;jt++){
        for(Index var=0;var<nY;var++){
          x[(Time+1)*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];
        }
        for (Index cup=0;cup<nU;cup++){
          x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];
        }
        for (Index par=0;par<nP;par++){
          x[(2*Time+1)*(nY+nU)+par] = rand()*1.0/RAND_MAX*(bounds[nY+nU+par][1]-bounds[nY+nU+par][0])+bounds[nY+nU+par][0];
        }
      }
    }
    for(Index i=0;i<2*Time+1;i++) delete [] init[i];
    delete [] init;
  }
  else{
    for(Index i=0;i<Ntotal-nP;i++){
      x[i] = solution[i];
    }
    for (Index jt=0;jt<Time+1;jt++){
      for (Index cup=0;cup<nU;cup++){
        x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];
        if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];
      }
    }
    for(Index i=0;i<nP;i++){
      x[(2*Time+1)*(nU+nY)+i] = solution[Ntotal-nP+i];
    }
  }
  return true;
}

// returns the value of the objective function
bool NEUTRINOS_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == 4*Time+2);
  obj_value = 0;

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[i*(Time+1) + jt];
        Xvalp1[i] = x[i*(Time+1) + jt + 1];
        Xval2[i] = x[(Time+1)*(nY+nU) + i*(Time) + jt];
     } //end for loop


     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + jt];
        K11valp1[cup] = x[nY*(Time+1) + cup*(Time+1) + jt + 1];
        K11val2[cup] = x[(Time+1)*(nY+nU) + Time*nY + cup*(Time) + jt];
     } //end for loop

     Xdval[0] = dd[2*jt];
     Xdval2[0] = dd[2*jt+1];
     Xdvalp1[0] = dd[2*jt+2];
     Xdval[1] = BC[2*jt];
     Xdval2[1] = BC[2*jt+1];
     Xdvalp1[1] = BC[2*jt+2];
     Xdval[2] = tt[2*jt];
     Xdval2[2] = tt[2*jt+1];
     Xdvalp1[2] = tt[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop


    obj_value += Xdval[1]*pow(-Xdval[0] + Xval[0], 2) + Xdval[1]*pow(Xval[1] - Xdval[2], 2) + Xdval2[1]*pow(-Xdval2[0] + Xval2[0], 2) + Xdval2[1]*pow(Xval2[1] - Xdval2[2], 2); 

    obj_value += bounds[0][3]*(pow(0.16666666666666699*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) + Xval[0] - Xvalp1[0], 2) + pow(0.125*hstep*(0.040000000000000001*pow(Xval[1], 3) - 0.040000000000000001*pow(Xvalp1[1], 3)) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0], 2)); 

    obj_value += bounds[1][3]*(pow(1.0*hstep + Xval[1] - Xvalp1[1], 2) + pow(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1], 2)); 

  } //end for loop

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[Time + i*(Time+1)];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop


     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[nY*(Time+1) + Time + cup*(Time+1)];
        K11valp1[cup] = 0;
        K11val2[cup] = 0;
     } //end for loop

     Xdval[0] = dd[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Xdval[1] = BC[2*Time];
     Xdval2[1] = 0;
     Xdvalp1[1] = 0;
     Xdval[2] = tt[2*Time];
     Xdval2[2] = 0;
     Xdvalp1[2] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop


  obj_value += Xdval[1]*pow(-Xdval[0] + Xval[0], 2) + Xdval[1]*pow(Xval[1] - Xdval[2], 2) + Xdval2[1]*pow(-Xdval2[0] + Xval2[0], 2) + Xdval2[1]*pow(Xval2[1] - Xdval2[2], 2);

  obj_value = obj_value/(2*Time+1);

  return true;
}


// return the gradient of the objective function grad_{x} f(x)
bool NEUTRINOS_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == 4*Time+2);

  for(Index i=0;i<n;i++) {
     grad_f[i] = 0;
  }

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[jt + i*(Time+1)];
        Xvalp1[i] = x[jt + i*(Time+1) + 1];
        Xval2[i] = x[(Time+1)*(nY+nU) + i*(Time)+ jt];
     } //end for loop

     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + jt];
        K11valp1[cup] = x[nY*(Time+1) + cup*(Time+1) + jt + 1];
        K11val2[cup] = x[(Time+1)*(nY+nU) + (nY)*Time + cup*(Time) +  jt];
     } //end for loop

     Xdval[0] = dd[2*jt];
     Xdval2[0] = dd[2*jt+1];
     Xdvalp1[0] = dd[2*jt+2];
     Xdval[1] = BC[2*jt];
     Xdval2[1] = BC[2*jt+1];
     Xdvalp1[1] = BC[2*jt+2];
     Xdval[2] = tt[2*jt];
     Xdval2[2] = tt[2*jt+1];
     Xdvalp1[2] = tt[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop

    grad_f[jt+0*(Time+1)] += (Xdval[1]*(-2*Xdval[0] + 2*Xval[0]))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += (Xdval[1]*(2*Xval[1] - 2*Xdval[2]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 0*Time + jt] += (Xdval2[1]*(-2*Xdval2[0] + 2*Xval2[0]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 1*Time + jt] += (Xdval2[1]*(2*Xval2[1] - 2*Xdval2[2]))/(2*Time+1);

    grad_f[jt+0*(Time+1)] += bounds[0][3]*(0.125*hstep*(0.040000000000000001*pow(Xval[1], 3) - 0.040000000000000001*pow(Xvalp1[1], 3)) + 0.33333333333333298*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) + 2.5*Xval[0] - 1.0*Xval2[0] - 1.5*Xvalp1[0])/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[0][3]*(0.040000000000000001*hstep*pow(Xval[1], 2)*(0.16666666666666699*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) + Xval[0] - Xvalp1[0]) + 0.029999999999999999*hstep*pow(Xval[1], 2)*(0.125*hstep*(0.040000000000000001*pow(Xval[1], 3) - 0.040000000000000001*pow(Xvalp1[1], 3)) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[0][3]*(0.125*hstep*(0.040000000000000001*pow(Xval[1], 3) - 0.040000000000000001*pow(Xvalp1[1], 3)) - 0.33333333333333298*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) - 1.5*Xval[0] - 1.0*Xval2[0] + 2.5*Xvalp1[0])/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[0][3]*(0.040000000000000001*hstep*pow(Xvalp1[1], 2)*(0.16666666666666699*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) + Xval[0] - Xvalp1[0]) - 0.029999999999999999*hstep*pow(Xvalp1[1], 2)*(0.125*hstep*(0.040000000000000001*pow(Xval[1], 3) - 0.040000000000000001*pow(Xvalp1[1], 3)) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 0*Time + jt] += bounds[0][3]*(-0.25*hstep*(0.040000000000000001*pow(Xval[1], 3) - 0.040000000000000001*pow(Xvalp1[1], 3)) - 1.0*Xval[0] + 2*Xval2[0] - 1.0*Xvalp1[0])/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 1*Time + jt] += bounds[0][3]*(0.16*hstep*pow(Xval2[1], 2)*(0.16666666666666699*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[1][3]*(2.0*hstep + 2.5*Xval[1] - 1.0*Xval2[1] - 1.5*Xvalp1[1])/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[1][3]*(-2.0*hstep - 1.5*Xval[1] - 1.0*Xval2[1] + 2.5*Xvalp1[1])/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 1*Time + jt] += bounds[1][3]*(-1.0*Xval[1] + 2*Xval2[1] - 1.0*Xvalp1[1])/(2*Time+1);

  } //end for loop

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[i*(Time+1) + Time];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop

     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + Time];
        K11valp1[cup] = 0;
        K11val2[cup] = 0;
     } //end for loop

     Xdval[0] = dd[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Xdval[1] = BC[2*Time];
     Xdval2[1] = 0;
     Xdvalp1[1] = 0;
     Xdval[2] = tt[2*Time];
     Xdval2[2] = 0;
     Xdvalp1[2] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop


    grad_f[Time+0*(Time+1)] += (Xdval[1]*(-2*Xdval[0] + 2*Xval[0]))/(2*Time+1);
    grad_f[Time+1*(Time+1)] += (Xdval[1]*(2*Xval[1] - 2*Xdval[2]))/(2*Time+1);

  return true;
}


// return the value of the constraints: g(x)
bool NEUTRINOS_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == 4*Time+2);
  assert(m == 0);

  return true;
}


// return the structure or values of the jacobian
bool NEUTRINOS_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index* jCol,
                            Number* values)
{

return true;
}


// return the structure or values of the hessian
bool NEUTRINOS_NLP::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{

if (values == NULL) {
   // return the structure.  This is a symmetric matrix, fill in the lower left
   // triangle only.

   // Each non-one Hessian element has its own explicit loop
   // since each element needs a different number of matrix elements

     for(Index jt=0;jt<Time+1;jt++) {
     iRow[0*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
     jCol[0*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[1*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt+1;
     jCol[1*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[1*(Time+1)+1*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[1*(Time+1)+1*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[2*(Time+1)+1*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[2*(Time+1)+1*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[2*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[2*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt+1;
     jCol[3*(Time+1)+2*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+3*(Time)+0+jt] = (Time+1)*1+jt+1;
     jCol[3*(Time+1)+3*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+4*(Time)+0+jt] = (Time+1)*2+Time*0+jt;
     jCol[3*(Time+1)+4*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+5*(Time)+0+jt] = (Time+1)*2+Time*0+jt;
     jCol[3*(Time+1)+5*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+6*(Time)+0+jt] = (Time+1)*2+Time*0+jt;
     jCol[3*(Time+1)+6*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+7*(Time)+0+jt] = (Time+1)*2+Time*0+jt;
     jCol[3*(Time+1)+7*(Time)+0+jt] = (Time+1)*1+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+8*(Time)+0+jt] = (Time+1)*2+Time*0+jt;
     jCol[3*(Time+1)+8*(Time)+0+jt] = (Time+1)*2+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+9*(Time)+0+jt] = (Time+1)*2+Time*1+jt;
     jCol[3*(Time+1)+9*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+10*(Time)+0+jt] = (Time+1)*2+Time*1+jt;
     jCol[3*(Time+1)+10*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+11*(Time)+0+jt] = (Time+1)*2+Time*1+jt;
     jCol[3*(Time+1)+11*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+12*(Time)+0+jt] = (Time+1)*2+Time*1+jt;
     jCol[3*(Time+1)+12*(Time)+0+jt] = (Time+1)*1+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+13*(Time)+0+jt] = (Time+1)*2+Time*1+jt;
     jCol[3*(Time+1)+13*(Time)+0+jt] = (Time+1)*2+Time*1+jt;
   }
}
else {
  // return the values.  This is a symmetric matrix, fill the lower left
  // triangle only
  // initialize the values array
  // Point to the initial starting spot for the Hessian elements

  for(Index jt=0;jt<3*(Time+1)+14*Time+0;jt++) values[jt] = 0.; // Initialize matrix

   // fill the objective portion

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[jt + i*(Time+1)];
        Xvalp1[i] = x[jt + i*(Time+1) + 1];
        Xval2[i] = x[(Time+1)*(nY+nU) + jt + i*(Time)];
     } //end for loop

     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[jt + nY*(Time+1) + cup*(Time+1)];
        K11valp1[cup] = x[jt + nY*(Time+1) + cup*(Time+1) + 1];
        K11val2[cup] = x[(Time+1)*(nY+nU) + (nY+cup)*Time + jt];
     } //end for loop

     Xdval[0] = dd[2*jt];
     Xdval2[0] = dd[2*jt+1];
     Xdvalp1[0] = dd[2*jt+2];
     Xdval[1] = BC[2*jt];
     Xdval2[1] = BC[2*jt+1];
     Xdvalp1[1] = BC[2*jt+2];
     Xdval[2] = tt[2*jt];
     Xdval2[2] = tt[2*jt+1];
     Xdvalp1[2] = tt[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*1*(2*Xdval[1])/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*1*(2*Xdval[1])/(2*Time+1);
   values[3*(Time+1)+8*(Time)+0+jt] += obj_factor*1*(2*Xdval2[1])/(2*Time+1);
   values[3*(Time+1)+13*(Time)+0+jt] += obj_factor*1*(2*Xdval2[1])/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[0][3]*(2.5)/(2*Time+1);
   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.5)/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(2.5)/(2*Time+1);
   values[1*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.055*hstep*pow(Xval[1], 2))/(2*Time+1);
   values[2*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.025000000000000001*hstep*pow(Xval[1], 2))/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.00125*pow(hstep, 2)*pow(Xval[1], 4) + 0.080000000000000002*hstep*Xval[1]*(0.16666666666666699*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) + Xval[0] - Xvalp1[0]) + 0.059999999999999998*hstep*Xval[1]*(0.125*hstep*(0.040000000000000001*pow(Xval[1], 3) - 0.040000000000000001*pow(Xvalp1[1], 3)) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[3*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.025000000000000001*hstep*pow(Xvalp1[1], 2))/(2*Time+1);
   values[1*(Time+1)+1*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.055*hstep*pow(Xvalp1[1], 2))/(2*Time+1);
   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.00035*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xvalp1[1], 2))/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.00125*pow(hstep, 2)*pow(Xvalp1[1], 4) + 0.080000000000000002*hstep*Xvalp1[1]*(0.16666666666666699*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) + Xval[0] - Xvalp1[0]) - 0.059999999999999998*hstep*Xvalp1[1]*(0.125*hstep*(0.040000000000000001*pow(Xval[1], 3) - 0.040000000000000001*pow(Xvalp1[1], 3)) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[3*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.0)/(2*Time+1);
   values[3*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.0)/(2*Time+1);
   values[3*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.029999999999999999*hstep*pow(Xval[1], 2))/(2*Time+1);
   values[3*(Time+1)+7*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.029999999999999999*hstep*pow(Xvalp1[1], 2))/(2*Time+1);
   values[3*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[0][3]*(2)/(2*Time+1);
   values[3*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.16*hstep*pow(Xval2[1], 2))/(2*Time+1);
   values[3*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.16*hstep*pow(Xval2[1], 2))/(2*Time+1);
   values[3*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0032000000000000002*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xval2[1], 2))/(2*Time+1);
   values[3*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0032000000000000002*pow(hstep, 2)*pow(Xval2[1], 2)*pow(Xvalp1[1], 2))/(2*Time+1);
   values[3*(Time+1)+13*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.012800000000000001*pow(hstep, 2)*pow(Xval2[1], 4) + 0.32000000000000001*hstep*Xval2[1]*(0.16666666666666699*hstep*(0.040000000000000001*pow(Xval[1], 3) + 0.16*pow(Xval2[1], 3) + 0.040000000000000001*pow(Xvalp1[1], 3)) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[1][3]*(2.5)/(2*Time+1);
   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.5)/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(2.5)/(2*Time+1);
   values[3*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.0)/(2*Time+1);
   values[3*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.0)/(2*Time+1);
   values[3*(Time+1)+13*(Time)+0+jt] += obj_factor*bounds[1][3]*(2)/(2*Time+1);
   } // end for loop 

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[Time + i*(Time+1)];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop

     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[Time + nY*(Time+1) + cup*(Time+1)];
        K11valp1[cup] = 0;
        K11val2[cup] = 0;
     } //end for loop

     Xdval[0] = dd[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Xdval[1] = BC[2*Time];
     Xdval2[1] = 0;
     Xdvalp1[1] = 0;
     Xdval[2] = tt[2*Time];
     Xdval2[2] = 0;
     Xdvalp1[2] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop


   values[0*(Time+1)+0*(Time)+0+ Time] += obj_factor*(2*Xdval[1])/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+ Time] += obj_factor*(2*Xdval[1])/(2*Time+1);
  } // end else 

   return true;
}


void NEUTRINOS_NLP::finalize_solution(SolverReturn status,
                        Index n, const Number* x, const Number* z_L, const Number* z_U,
                        Index m, const Number* g, const Number* lambda,
                        Number obj_value,
                        const IpoptData* ip_data,
                        IpoptCalculatedQuantities* ip_cq)
{
  FILE *OUTPUT1;
  char filename[20];
  if ((file_fmt < 2) || (modid == 0)){
    sprintf(filename,"D%d_M%d_IC%d.dat", nY,nM,taskid);
  }
  else if ((file_fmt >= 2) && (modid != 0)){
    sprintf(filename,"D%d_M%d_PATH%d_IC%d.dat", nY,nM,pathid,taskid);
  }
  if(beta==0){
    OUTPUT1 = fopen (filename,"w");
  }
  else{
    if ((output_fmt == 1) || (output_fmt < 0)){
      OUTPUT1 = fopen (filename,"w");
    }
    else 
      OUTPUT1 = fopen (filename,"a");
  }
  // Write solution for annealing
  for(Index i=0;i<Ntotal-nP;i++){
    solution[i] = x[i];
  }
  for(Index i=0;i<nP;i++){
    solution[Ntotal-nP+i] = x[(2*Time+1)*(nU+nY)+i];
  }

  // Write to file
  if (output_fmt != 1){
    fprintf(OUTPUT1, "%d %d %e ",beta, (status == SUCCESS), obj_value);
  }
  for (Index i=0;i<Time;i++) {
    for (Index j=0;j<nY;j++) {
      fprintf(OUTPUT1,"%e ", x[j*(Time+1)+i]);
    }
    if (abs(output_fmt) == 2){
      for (Index cup=0;cup<nU;cup++) {
        fprintf(OUTPUT1,"%e ", x[nY*(Time+1)+ cup*(Time+1)+ i]);
      }
    }
    for (Index j=0;j<nY;j++) {
      fprintf(OUTPUT1,"%e ", x[(nY+nU)*(Time+1) + j*Time + i]);
    }
    if (abs(output_fmt) == 2){
      for (Index cup=0;cup<nU;cup++) {
        fprintf(OUTPUT1,"%e ", x[(nY+nU)*(Time+1) + nY*Time + cup*Time + i]);
      }
    }
  }
  for (Index j=0;j<nY;j++) {
     fprintf(OUTPUT1,"%e ", x[j*(Time+1) + Time]);
  }
  if (abs(output_fmt) == 2){
    for (Index cup=0;cup<nU;cup++) {
       fprintf(OUTPUT1,"%e ", x[nY*(Time+1) + cup*(Time+1) + Time]);
    }
  }
  for (Index j=0;j<nP;j++) {
     fprintf(OUTPUT1,"%e ", x[(2*Time+1)*(nY+nU)+j]);
  }
  fprintf(OUTPUT1,"\n");
  fclose (OUTPUT1);
  

}
