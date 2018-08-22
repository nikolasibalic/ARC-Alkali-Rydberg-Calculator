
#include <Python.h>
// http://docs.scipy.org/doc/numpy/reference/c-api.deprecations.html
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double innerLimit;
double outerLimit;
double step;
double init1;
double init2;

int l;
double s,j;
double stateEnergy;
double alphaC;
double alpha;
int Z;
double a1,a2,a3,a4,rc; // depends on l - determined by Python in advance
double mu;

//#define DEBUG_OUTPUT

static PyObject *NumerovWavefunction(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
   {"NumerovWavefunction", NumerovWavefunction, METH_VARARGS,
    "Numerov wavefunction"},
  {NULL, NULL, 0, NULL}};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT, "arc_c_extensions",
  "C extensions of ARC (Numerov integration)", -1, module_methods, NULL, NULL, NULL, NULL, };

PyMODINIT_FUNC PyInit_arc_c_extensions(void) {
  // import Numpy API
  import_array();

  return PyModule_Create(&moduledef);
}
#else
PyMODINIT_FUNC initarc_c_extensions(void) {
  if (!(Py_InitModule3("arc_c_extensions", module_methods,
                       "C extensions of ARC (Numerov integration)"))) return;
  //  import Numpy API
  import_array();
}
#endif


// =========== variable definition ===========
int divergencePoint;
int totalLength;
int br ;
double* sol;
double x,step2,maxValue,checkPoint,fromLastMax,r;
int i;
npy_intp dims[2];
PyObject* narray;

double commonTerm1=0;
double commonTerm2=0;

// =========== Numerov integration implementation ===========

__inline double EffectiveCharge(double r){
	// returns effective charge of the core
	return 1.0+((double)Z-1.)*exp(-a1*r)-r*(a3+a4*r)*exp(-a2*r);
}

__inline  double CorePotential(double r){
    // l dependent core potential (angular momentum of e) at radius r
    // returns Core potential
    return -EffectiveCharge(r)/r-alphaC/(2.*pow(r,4))*(1.-exp(-pow(r/rc,6)));
}


__inline double Potenital(double r){
	// l<4
    return CorePotential(r)+pow(alpha,2)/(2.0*pow(r,3))*commonTerm1;

}

__inline double Potenital2(double r){
	// l>=4
    // act as if it is a Hydrogen atom, include spin-orbit coupling
    return -1./r+pow(alpha,2)/(2.0*pow(r,3))*commonTerm1;
}


__inline double kfun(double x){
	// with potential for l<4
  r = x*x;   // x = sqrt(r)
	return -3./(4.*r)+4*r*( 2.*mu*(stateEnergy-Potenital(r))-commonTerm2/pow(r,2) );
}

__inline double kfun2(double x){
	// with potential for l>4
  r = x*x;  // x = sqrt(r)
	return -3./(4.*r)+4*r*( 2.*mu*(stateEnergy-Potenital2(r))-commonTerm2/pow(r,2) );
}

static PyObject *NumerovWavefunction(PyObject *self, PyObject *args) {
	// Numerov arguments: innerLimit,outerLimit,kfun,step,init1,init2
	double innerLimit,outerLimit,step,init1,init2;


    if (!(PyArg_ParseTuple(args, "dddddidddddidddddd", &innerLimit, &outerLimit, &step,
      &init1, &init2,
      &l, &s, &j, &stateEnergy, &alphaC,  &alpha,
      &Z, &a1, &a2, &a3, &a4, &rc, &mu))) return NULL;


#ifdef DEBUG_OUTPUT
		printf("innerLimit\t=\t%.3f\nouterLimit\t=\t%.3f\nstep\t=\t%.3f\ninit1\t=\t%.3f\ninit2\t=\t%.3f\n",innerLimit,outerLimit,step,init1,init2);
		printf("l\t=\t%i\ns\t=\t%.1f\nj\t=\t%.1f\n",l,s,j);
		printf("stateEnergy\t=\t%.7f\nalphaC\t=\t%.3f\nalpha\t=\t%.3f\nZ\t=\t%i\n",stateEnergy,alphaC,alpha,Z);
		printf("a1\t\t%.4f\na2\t\t%.4f\na3\t\t%.4f\na4\t\t%.4f\nrc\t\t%.4f\n",a1,a2,a3,a4,rc);
    printf("mu\t\t%.4f",mu);
#endif

	// let's speed up calculation by calculating some common terms beforehand
	commonTerm1 = (j*(j+1.0)-((double)l)*(l+1.0)-s*(s+1.))/2.0;
	commonTerm2 = ((double)l)*(l+1.);

	totalLength =  (int)((sqrt(outerLimit)-sqrt(innerLimit))/step);

#ifdef DEBUG_OUTPUT
	printf("Index = %i\n",totalLength);
	printf("Index should be about = %.2f\n",(sqrt(outerLimit)-sqrt(innerLimit)/step));
#endif

	br = totalLength;
	sol = (double*) malloc(2*br*sizeof(double));

	if (!sol){
  #ifdef DEBUG_OUTPUT
		printf("Memory allocaiton failed! Aborting.");
  #endif
		return NULL;
	}

    // for l<4

	if (l<4){

		br = br-1;
    x = sqrt(innerLimit)+step*(totalLength-1);
	  step2 = step*step;
	  sol[br] = (2*(1-5.0/12.0*step2*kfun(x))*init1-(1+1/12.0*step2*kfun(x+step))*init2)/(1+1/12.0*step2*kfun(x-step));
    sol[br+totalLength]=x;
		x = x-step;
		br = br-1;

		sol[br] = (2*(1-5.0/12.0*step2*kfun(x))*sol[br+1]-(1+1/12.0*step2*kfun(x+step))*init1)/(1+1/12.0*step2*kfun(x-step));
    sol[br+totalLength]=x;

		maxValue = 0;

	  checkPoint = 0;
	  fromLastMax = 0;

	  while (br>checkPoint){
	      br = br-1;
	      x = x-step;
	      sol[br] = (2*(1-5.0/12.0*step2*kfun(x))*sol[br+1]-(1+1/12.0*step2*kfun(x+step))*sol[br+2])/(1+1/12.0*step2*kfun(x-step));
        sol[br+totalLength]=x;
	      if (fabs(sol[br]*sqrt(x))>maxValue){
	          maxValue = fabs(sol[br]*sqrt(x));
        }
        else{
	          fromLastMax += 1;
            if (fromLastMax>50){
	              checkPoint = br;
	          }
	      }
	  }

	  divergencePoint = 0;
	  while ((br>0)&&(divergencePoint == 0)){
	      br = br-1;
	      x = x-step;
	      sol[br] = (2*(1-5.0/12.0*step2*kfun(x))*sol[br+1]-(1+1/12.0*step2*kfun(x+step))*sol[br+2])/(1+1/12.0*step2*kfun(x-step));
        sol[br+totalLength]=x;

	      if ((divergencePoint==0)&&(fabs(sol[br]*sqrt(x))>maxValue)){
	          divergencePoint = br;
	          while ((fabs(sol[divergencePoint])>fabs(sol[divergencePoint+1])) && (divergencePoint<checkPoint)){
	              divergencePoint +=1;
	          }
	          if (divergencePoint>checkPoint){
#ifdef DEBUG_OUTPUT
	                printf("ERROR: Numerov error\n");
#endif
	                return NULL;
	          }
	      }
	  }


	} // end of if l<4
	else{ //l>=4

		  br = br-1;
	    x = sqrt(innerLimit)+step*(totalLength-1);
	    step2 = step*step;
	    sol[br] = (2*(1-5.0/12.0*step2*kfun2(x))*init1-(1+1/12.0*step2*kfun2(x+step))*init2)/(1+1/12.0*step2*kfun2(x-step));
      sol[br+totalLength]=x;
		  x = x-step;
		  br = br-1;

		  sol[br] = (2*(1-5.0/12.0*step2*kfun2(x))*sol[br+1]-(1+1/12.0*step2*kfun2(x+step))*init1)/(1+1/12.0*step2*kfun2(x-step));
      sol[br+totalLength]=x;

		  maxValue = 0;

	    checkPoint = 0;
	    fromLastMax = 0;

	    while (br>checkPoint){
	        br = br-1;
	        x = x-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun2(x))*sol[br+1]-(1+1/12.0*step2*kfun2(x+step))*sol[br+2])/(1+1/12.0*step2*kfun2(x-step));
          sol[br+totalLength]=x;

    	    if (fabs(sol[br]*sqrt(x))>maxValue){
    	        maxValue = fabs(sol[br]*sqrt(x));
    	    }
    	    else{
    	         fromLastMax += 1;
    	         if (fromLastMax>50){
    	             checkPoint = br;
    	         }
	        }
	    }

	    divergencePoint = 0;
	    while ((br>0)&&(divergencePoint == 0)){
	        br = br-1;
	        x = x-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun2(x))*sol[br+1]-(1+1/12.0*step2*kfun2(x+step))*sol[br+2])/(1+1/12.0*step2*kfun2(x-step));
          sol[br+totalLength]=x;

	        if ((divergencePoint==0)&&(fabs(sol[br]*sqrt(x))>maxValue)){
	            divergencePoint = br;
	            while ((fabs(sol[divergencePoint])>fabs(sol[divergencePoint+1])) && (divergencePoint<checkPoint)){
	                divergencePoint +=1;
	            }
	            if (divergencePoint>checkPoint){
#ifdef DEBUG_OUTPUT
	                printf("ERROR: Numerov error\n");
#endif
	                return NULL;
	            }
	        }
	    }

	}

  // RETURN RESULT - but set to zero divergent part (to prevent integration there)
  for (i =0; i<divergencePoint; i++) sol[i] = 0;
  // same for radial part
  for (i = divergencePoint; i >= 0 ; i--) sol[i+totalLength] = sol[i+totalLength+1]-step;

  // convert sol that is at the moment R(r)*r^{3/4} into R(r)*r
  for (i=0; i<totalLength; i++)  sol[i]=sol[i]*sqrt(sol[i+totalLength]);
  // convert coordinates from sqrt(r) into r
  for (i=totalLength; i<2*totalLength; i++)  sol[i]=sol[i]*sol[i];

  // return the array as a numpy array (numpy will free it later)
  dims[0] = totalLength;
  dims[1] = totalLength;
  narray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sol);
  //free(sol); # freeing of solution array should be done from Numpy
  // this is the critical line - tell numpy it has to free the data
  PyArray_ENABLEFLAGS((PyArrayObject*)narray, NPY_ARRAY_OWNDATA);
  return narray;

  return 0;
}
