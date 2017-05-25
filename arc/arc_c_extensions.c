
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

//#define DEBUG_OUTPUT

static PyObject *NumerovWavefunction(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
   {"NumerovWavefunction", NumerovWavefunction, METH_VARARGS,
    "Numerov wavefunction"},
  {NULL, NULL, 0, NULL}};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT, "arc_c_extensions", 
  "C extensions of ARC (Numerov integration)", -1, module_methods, };

PyMODINIT_FUNC PyInit_arc_c_extensions(void) {
  return PyModule_Create(&moduledef);

  // something to do with numpy
  import_array();
}
#else
PyMODINIT_FUNC initarc_c_extensions(void) {
  if (!(Py_InitModule3("arc_c_extensions", module_methods,
                       "C extensions of ARC (Numerov integration)"))) return;
  // something to do with numpy
  import_array();
}
#endif





// numerov

inline double EffectiveCharge(double r){
	// returns effective charge of the core
	return 1.0+((double)Z-1.)*exp(-a1*r)-r*(a3+a4*r)*exp(-a2*r);
}

inline double CorePotential(double r){
    // l dependent core potential (angular momentum of e) at radius r
    // returns Core potential
    return -EffectiveCharge(r)/r-alphaC/(2.*pow(r,4))*(1.-exp(-pow(r/rc,6)));
}

double commonTerm1;

inline double Potenital(double r){
	// l<4
    return CorePotential(r)+pow(alpha,2)/(2.0*pow(r,3))*commonTerm1;

}

inline double Potenital2(double r){
	// l>=4
    // act as if it is a Hydrogen atom
    return -1./r;
}

double commonTerm2;

inline double kfun(double r){
	// with potential for l<4
	return 2.*(stateEnergy-Potenital(r))-commonTerm2/pow(r,2);
}

inline double kfun2(double r){
	// with potential for l>4
	return 2.*(stateEnergy-Potenital2(r))-commonTerm2/pow(r,2);
}

static PyObject *NumerovWavefunction(PyObject *self, PyObject *args) {
	// Numerov arguments: innerLimit,outerLimit,kfun,step,init1,init2
	double innerLimit,outerLimit,step,init1,init2;


    if (!(PyArg_ParseTuple(args, "dddddidddddiddddd", &innerLimit, &outerLimit, &step,
      &init1, &init2,
      &l, &s, &j, &stateEnergy, &alphaC,  &alpha,
      &Z, &a1, &a2, &a3, &a4, &rc))) return NULL;


#ifdef DEBUG_OUTPUT
		printf("innerLimit\t=\t%.3f\nouterLimit\t=\t%.3f\nstep\t=\t%.3f\ninit1\t=\t%.3f\ninit2\t=\t%.3f\n",innerLimit,outerLimit,step,init1,init2);
		printf("l\t=\t%i\ns\t=\t%.1f\nj\t=\t%.1f\n",l,s,j);
		printf("stateEnergy\t=\t%.7f\nalphaC\t=\t%.3f\nalpha\t=\t%.3f\nZ\t=\t%i\n",stateEnergy,alphaC,alpha,Z);
		printf("a1\t\t%.4f\na2\t\t%.4f\na3\t\t%.4f\na4\t\t%.4f\nrc\t\t%.4f\n",a1,a2,a3,a4,rc);
#endif


	// let's speed up calculation by calculating some common terms beforehand
	commonTerm1 = (j*(j+1.0)-((double)l)*(l+1.0)-s*(s+1.))/2.0;
	commonTerm2 = ((double)l)*(l+1.);

	int divergencePoint;

	int totalLength =  (int)((outerLimit-innerLimit)/step);

#ifdef DEBUG_OUTPUT
	printf("Index = %i\n",totalLength);
	printf("Index should be about = %.2f\n",((outerLimit-innerLimit)/step));
#endif

	int br = totalLength;
	double* sol = (double*) malloc(br*sizeof(double));

	if (!sol){
  #ifdef DEBUG_OUTPUT
		printf("Memory allocaiton failed! Aborting.");
  #endif
		return NULL;
	}

    // for l<4

	if (l<4){

		br = br-1;
	    double r = outerLimit;
	    double step2 = step*step;
	    sol[br] = (2*(1-5.0/12.0*step2*kfun(r))*init1-(1+1/12.0*step2*kfun(r+step))*init2)/(1+1/12.0*step2*kfun(r-step));

		r = r-step;
		br = br-1;

		sol[br] = (2*(1-5.0/12.0*step2*kfun(r))*sol[br+1]-(1+1/12.0*step2*kfun(r+step))*init1)/(1+1/12.0*step2*kfun(r-step));

		double maxValue = 0;

	    double checkPoint = 0;
	    double fromLastMax = 0;

	    while (br>checkPoint){
	        br = br-1;
	        r = r-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun(r))*sol[br+1]-(1+1/12.0*step2*kfun(r+step))*sol[br+2])/(1+1/12.0*step2*kfun(r-step));
	        if (fabs(sol[br])>maxValue){
	            maxValue = fabs(sol[br]);
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
	        r = r-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun(r))*sol[br+1]-(1+1/12.0*step2*kfun(r+step))*sol[br+2])/(1+1/12.0*step2*kfun(r-step));
	        if ((divergencePoint==0)&&(fabs(sol[br])>maxValue)){
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
	    double r = outerLimit;
	    double step2 = step*step;
	    sol[br] = (2*(1-5.0/12.0*step2*kfun2(r))*init1-(1+1/12.0*step2*kfun2(r+step))*init2)/(1+1/12.0*step2*kfun2(r-step));

		r = r-step;
		br = br-1;

		sol[br] = (2*(1-5.0/12.0*step2*kfun2(r))*sol[br+1]-(1+1/12.0*step2*kfun2(r+step))*init1)/(1+1/12.0*step2*kfun2(r-step));

		double maxValue = 0;

	    double checkPoint = 0;
	    double fromLastMax = 0;

	    while (br>checkPoint){
	        br = br-1;
	        r = r-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun2(r))*sol[br+1]-(1+1/12.0*step2*kfun2(r+step))*sol[br+2])/(1+1/12.0*step2*kfun2(r-step));
	        if (fabs(sol[br])>maxValue){
	            maxValue = fabs(sol[br]);
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
	        r = r-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun2(r))*sol[br+1]-(1+1/12.0*step2*kfun2(r+step))*sol[br+2])/(1+1/12.0*step2*kfun2(r-step));
	        if ((divergencePoint==0)&&(fabs(sol[br])>maxValue)){
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
  for (int i =0; i<divergencePoint; i++) sol[i] = 0;

  // return the array as a numpy array (numpy will free it later)
  npy_intp dims[1] = {totalLength};
  PyObject *narray = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, sol);
  //free(sol); # freeing of solution array should be done from Numpy
  // this is the critical line - tell numpy it has to free the data
  PyArray_ENABLEFLAGS((PyArrayObject*)narray, NPY_ARRAY_OWNDATA);
  return narray;

  return 0;
}
