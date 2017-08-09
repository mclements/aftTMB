#include <TMB.hpp>      // Links in the TMB libraries

// typedef matrix<double> MatrixXd;
// typedef vector<double> VectorXd;
// // conversions
// template<class Type>
// matrix<Type> asMatrix(MatrixXd a) {
//   matrix<Type> b(a.rows(), a.cols());
//   for (int i=0; i<a.rows(); i++)
// 	for (int j=0; j<a.cols(); j++)
// 	  b(i,j) = Type(a(i,j));
//   return b;
// }
// template<class Type>
// MatrixXd asMatrixXd(matrix<Type> a) {
//   MatrixXd b(a.rows(), a.cols());
//   for (int i=0; i<a.rows(); i++)
// 	for (int j=0; j<a.cols(); j++)
// 	  b(i,j) = asDouble(a(i,j));
//   return b;
// }
// template<class Type>
// vector<Type> asVector(VectorXd a) {
//   vector<Type> b(a.size());
//   for (int i=0; i<a.size(); i++)
// 	b(i) = Type(a(i));
//   return b;
// }
// template<class Type>
// VectorXd asVectorXd(vector<Type> a) {
//   VectorXd b(a.size());
//   for (int i=0; i<a.size(); i++)
// 	b(i) = asDouble(a(i));
//   return b;
// }
// template<class Type>
// matrix<Type> head(MatrixXd a, int n = 6) {
//   matrix<Type> b(n, a.cols());
//   for (int i=0; i<n; i++)
// 	for (int j=0; j<a.cols(); j++)
// 	  b(i,j) = a(i,j);
//   return b;
// }
// template<class Type>
// vector<Type> head(vector<Type> a, int n = 6) {
//   vector<Type> b(n);
//   for (int i=0; i<n; i++)
// 	b(i) = a(i);
//   return b;
// }
// dot products
template<class Type>
Type dot(vector<Type> a, vector<Type> b) {
  return (a.array() * b.array()).sum();
}
template<class Type>
Type dot3(vector<Type> a, vector<Type> b, vector<Type> c) {
  return (a.array() * b.array() * c.array()).sum();
}
// // printing vectors and matrices
// template<class Type>
// void Rprint(vector<Type> a) {
//   for (int i=0; i<a.size(); i++) {
// 	Rprintf(i==0 ? "[" : ",");
// 	Rprintf("%f",asDouble(a(i)));
//   }
//   Rprintf("]\n");
// }
// template<class Type>
// void Rprint(matrix<Type> a) {
//   for (int i=0; i<a.rows(); i++) {
// 	Rprintf(i==0 ? "[" : " ");
// 	for (int j=0; j<a.cols(); j++) {
// 	  Rprintf(j==0 ? "[" : ",");
// 	  Rprintf("%f",asDouble(a(i,j)));
// 	}
// 	Rprintf("]\n");
//   }
//   Rprintf("]\n");
// }
// // Complete Q matrix from a QR decomposition
// matrix<double> qr_q(const matrix<double>& X, double tol = 1E-12) 
// {
//   // Initialize member data and allocate heap memory
//   int n=X.rows(), p=X.cols(), rank=0;
//   matrix<double> qr = X, y(n,n), q(n,n);
//   int* pivot=(int*)R_alloc(p,sizeof(int)); 
//   double* tau=(double*)R_alloc(p,sizeof(double)); 
//   double* work=(double*)R_alloc(p*2,sizeof(double));
//   for(int i=0;i<p;i++) 
// 	pivot[i]=i+1; 
//   y.setIdentity();
//   // LINPACK QR factorization via householder transformations
//   F77_CALL(dqrdc2)(qr.data(), &n, &n, &p, &tol, &rank, tau, pivot, work);
//   // Compute orthogonal factor Q
//   F77_CALL(dqrqy)(qr.data(), &n, &rank, tau, y.data(), &n, q.data());
//   return q;
// }

template<class Type>
int asInt(Type x) { return (int) asDouble(x); }
template<class Type>
class SplineBasis {
public:
  int order,			/* order of the spline */
    ordm1,			/* order - 1 (3 for cubic splines) */
    nknots,			/* number of knots */
    ncoef;			/* number of coefficients */
  Type curs,			/* current position in knots vector */
    boundary,			/* must have knots[(curs) <= x < knots(curs+1) */
    offset,			/* offset */
    branch;			/* branch indicator */
  /* except for the boundary case */
  vector<Type> ldel;  	/* differences from knots on the left */
  vector<Type> rdel;	/* differences from knots on the right */
  vector<Type> knots;	/* knot vector */
  vector<Type> a;		/* scratch array */
  SplineBasis(int order = 4) : order(order) {
    ordm1 = order - 1;
    rdel = vector<Type>(ordm1);
    ldel = vector<Type>(ordm1);
    a = vector<Type>(order);
    offset = 0;
  }
  SplineBasis(vector<Type> knots, int order = 4) : order(order), knots(knots) {
    ordm1 = order - 1;
    nknots = knots.size();
    ncoef = nknots - order;
    rdel = vector<Type>(ordm1);
    ldel = vector<Type>(ordm1);
    a = vector<Type>(order);
    offset = 0;
  }
  Type set_cursor(Type x)
  {
    Type knot, test_gt;
    int i;
    curs = Type(-1); /* Wall */
    boundary = Type(0);
    int lastLegit = nknots - order;
    test_gt = Type(0);
    for (i = 0; i < nknots; i++) {
      knot = knots(i);
      curs = CppAD::CondExpGe(knot, x, CppAD::CondExpEq(test_gt, Type(1), curs, Type(i)), curs);
      test_gt = CppAD::CondExpGt(knot,x,Type(1),test_gt);
    }
    knot = knots(lastLegit);
    boundary = CppAD::CondExpGt(curs, Type(nknots - order),
				CppAD::CondExpEq(x,knot,Type(1),Type(0)), Type(0));
    curs = CppAD::CondExpGt(curs, Type(nknots - order),
			    CppAD::CondExpEq(x,knot,Type(lastLegit),curs), curs);
    return curs;
  }
  void
  diff_table(Type x, int ndiff)
  {
    for (int i = 0; i < ndiff; i++) {
      rdel(i) = knots(asInt(curs) + i) - x;
      ldel(i) = x - knots(asInt(curs) - (i + 1));
    }
  }
  Type slow_evaluate(Type x, int i, int nder)
  {
    Type ti = curs, lpt;
    int  
      apt, rpt, inner,
      outer = ordm1;
    for(int j = 0; j < order; j++)
      a(j) = Type(0);
    a(i) = Type(1);
    //if (boundary && nder == ordm1) { /* value is arbitrary */
    //  return Type(0);
    //}
    while(nder--) {  // FIXME: divides by zero
      for(inner = outer, apt = 0, lpt = ti - outer; inner--; apt++, lpt+=1)
	a(apt) = Type(outer) * (a(apt + 1) - a(apt))/(knots(asInt(lpt) + outer) - knots(asInt(lpt)));
      outer--;
    }
    diff_table(x, outer);
    while(outer--)
      for(apt = 0, lpt = outer, rpt = 0, inner = outer + 1;
	  inner--; lpt-=1, rpt++, apt++)
	// FIXME: divides by zero
	a(apt) = (a(apt + 1) * ldel(asInt(lpt)) + a(apt) * rdel(rpt))/(rdel(rpt) + ldel(asInt(lpt)));
    return a(0);
  }
  /* fast evaluation of basis functions - NOT used in evali() */
  vector<Type> basis_funcs(Type x)
  {
    vector<Type> b(order);
    set_cursor(x);
    diff_table(x, ordm1);
    b(0) = Type(1);
    for (int j = 1; j <= ordm1; j++) {
      Type saved = Type(0);
      for (int r = 0; r < j; r++) { // do not divide by zero
	Type den = rdel(r) + ldel(j - 1 - r);
	Type term = CppAD::CondExpEq(den, Type(0), Type(0), b(r)/den);
	b(r) = CppAD::CondExpEq(den, Type(0),
				CppAD::CondExpEq(Type(r), Type(0),
						 CppAD::CondExpEq(rdel(r), Type(0), b(r), saved), saved), saved + rdel(r) * term);
	saved = CppAD::CondExpEq(den, Type(0), Type(0), ldel(j - 1 - r) * term);
      }
      b(j) = saved;
    }
    return b;
  }
  void set_options(Type x) {
    set_cursor(x);
    offset = curs - Type(order);
    // set branch=-1 if offset is out of range, otherwise branch=1
    branch = CppAD::CondExpLt(offset,Type(0),Type(-1),
			      CppAD::CondExpGt(offset,Type(nknots),Type(-1),Type(1)));
    // reset offset=0 if offset is out of range
    offset = CppAD::CondExpEq(branch,Type(-1),Type(0),offset); 
  }
  Type evali(Type x, int i, int ders=0) {
    Type val;
    val = CppAD::CondExpEq(branch,Type(-1),Type(0),
			   slow_evaluate(x, i, ders));
    return val;
  }
  matrix<Type> basis(vector<Type> x, int ders = 0) {
    matrix<Type> mat(x.size(), ncoef);
    mat.setZero();
    for (int i=0; i<x.size(); i++) {
      set_options(x(i));
      for  (int j=0; j<order; j++)
	mat(i,j+asInt(offset))=evali(x(i), j, ders);
    }
    return mat;
  }
};
template<class Type>
class bs : public SplineBasis<Type> {
public:
  vector<Type> boundary_knots, interior_knots;
  int intercept, df;
  bs(vector<Type> boundary_knots, vector<Type> interior_knots, int intercept = 0) :
    SplineBasis<Type>(4), boundary_knots(boundary_knots), interior_knots(interior_knots),
    intercept(intercept) {
    df = intercept + 3 + interior_knots.size();
    this->nknots = interior_knots.size()+8;
    this->ncoef = this->nknots - this->order;
    this->knots = vector<Type>(this->nknots);
    for(int i=0; i<4;i++) {
      this->knots(i)=boundary_knots(0);
      this->knots(this->nknots-i-1)=boundary_knots(1);
    }
    if (interior_knots.size() > 0) 
      for(int i=0; i<interior_knots.size();i++) 
	this->knots(i+4)=interior_knots(i);
  }
  vector<Type> eval(Type x, int ders=0) {
    vector<Type> xv(1);
    xv(0) = x;
    matrix<Type> mat = bs<Type>::basis(xv, ders);
    vector<Type> val(mat.cols());
    for (int i=0; i<mat.cols(); i++)
      val(i) = mat(0,i);
    return val;
  }
  matrix<Type> basis(vector<Type> x, int ders=0) {
    matrix<Type> mat(x.size(), this->df);
    mat.setZero();
    Type lower=boundary_knots(0), upper=boundary_knots(1);
    for (int i=0; i<x.size(); i++) {
      Type k_pivot = CppAD::CondExpLt(x(i),lower,Type(0.75)*lower+Type(0.25)*interior_knots(0), 
				      CppAD::CondExpGt(x(i),upper,
						       Type(0.75)*upper+Type(0.25)*interior_knots(interior_knots.size()-1),
						       x(i)));
      Type delta = x(i) - k_pivot; // not needed between the boundaries
      this->set_options(k_pivot);
      for  (int j=1-intercept; j < this->order; j++) {
	Type val;
	val = CppAD::CondExpLt(x(i), lower, 
			       SplineBasis<Type>::evali(k_pivot,j,0) +
			       SplineBasis<Type>::evali(k_pivot,j,1)*delta +
			       SplineBasis<Type>::evali(k_pivot,j,2)*delta*delta/Type(2) +
			       SplineBasis<Type>::evali(k_pivot,j,3)*delta*delta*delta/Type(6),
			       CppAD::CondExpGt(x(i), upper, 
						SplineBasis<Type>::evali(k_pivot,j,0) +
						SplineBasis<Type>::evali(k_pivot,j,1)*delta +
						SplineBasis<Type>::evali(k_pivot,j,2)*delta*delta/Type(2) +
						SplineBasis<Type>::evali(k_pivot,j,3)*delta*delta*delta/Type(6),
						SplineBasis<Type>::evali(k_pivot, j, ders)));
	mat(i,j+asInt(this->offset)-1+intercept)=val;
      }
    }
    return mat;
  }
};
#define ADIF(var,valt,valf) CppAD::CondExpEq(var,Type(1),valt,valf)
#define ADLT(var,val) CppAD::CondExpLt(var,val,Type(1),Type(0))
#define ADGT(var,val) CppAD::CondExpGt(var,val,Type(1),Type(0))
#define ADEQ(var,val) CppAD::CondExpEq(var,val,Type(1),Type(0))
#define ADGE(var,val) CppAD::CondExpGe(var,val,Type(1),Type(0))
#define ADLE(var,val) CppAD::CondExpLe(var,val,Type(1),Type(0))
#define ADNOT(var) CppAD::CondExpEq(var,Type(1),Type(0),Type(1))
#define ADOR(cond1,cond2,valt,valf) CppAD::CondExpEq(cond1,Type(1),valt,CppAD::CondExpEq(cond2,Type(1),valt,valf))
#define ADAND(cond1,cond2,valt,valf) CppAD::CondExpEq(cond1,Type(1),CppAD::CondExpEq(cond2,Type(1),valt,valf), valf)
template<class Type>
class ns : public bs<Type> {
public:
  vector<Type> tl0, tl1, tr0, tr1;
  matrix<Type> q_matrix;
  ns(vector<Type> boundary_knots, vector<Type> interior_knots, matrix<Type> _q_matrix, int intercept=0) :
    bs<Type>(boundary_knots, interior_knots, intercept), q_matrix(_q_matrix) {
    tl0 = q_matrix * bs<Type>::eval(boundary_knots(0), 0);
    tl1 = q_matrix * bs<Type>::eval(boundary_knots(0), 1);
    tr0 = q_matrix * bs<Type>::eval(boundary_knots(1), 0);
    tr1 = q_matrix * bs<Type>::eval(boundary_knots(1), 1);
  }
  matrix<Type> basis(vector<Type> x, int der=0) {
    matrix<Type> mat(x.size(), this->df-2);
    Type cl0, cl1, cr0, cr1, ci, k_pivot;
    for (int i=0; i<x.size(); i++) {
      ci = CppAD::CondExpLt(x(i),this->boundary_knots(0),Type(0),
			    CppAD::CondExpGt(x(i),this->boundary_knots(1), Type(0), Type(1)));
      if (der==0) {
	cl0 = CppAD::CondExpLt(x(i),this->boundary_knots(0),Type(1), Type(0));
	cl1 = CppAD::CondExpLt(x(i),this->boundary_knots(0),x(i)-this->boundary_knots(0),Type(0));
	cr0 = CppAD::CondExpGt(x(i),this->boundary_knots(1),Type(1), Type(0));
	cr1 = CppAD::CondExpGt(x(i),this->boundary_knots(1),x(i)-this->boundary_knots(1),Type(0));
      } else if (der==1) {
	cl0 = Type(0);
	cl1 = CppAD::CondExpLt(x(i),this->boundary_knots(0),Type(1),Type(0));
	cr0 = Type(0);
	cr1 = CppAD::CondExpGt(x(i),this->boundary_knots(1),Type(1),Type(0));
      } else {
	cl0 = Type(0);
	cl1 = Type(0);
	cr0 = Type(0);
	cr1 = Type(0);
      }
      k_pivot = CppAD::CondExpLt(x(i),this->boundary_knots(0),this->boundary_knots(0),CppAD::CondExpGt(x(i),this->boundary_knots(1),this->boundary_knots(1), x(i)));
      // Rprintf("x(%i)=%f,k_pivot=%f\n",i,asDouble(x(i)),k_pivot); 
      vector<Type> val = q_matrix * bs<Type>::eval(k_pivot,der);
      for  (int j=0; j<val.size(); j++)
	mat(i,j)=cl0*tl0(j)+cl1*tl1(j)+cr0*tr0(j)+cr1*tr1(j)+ci*val(j);
    }
    return mat;
  }
}; // class ns

template<class Type>
vector<Type> pmin(vector<Type> left, vector<Type> right) {
  vector<Type> val(left.size());
  for(int i=0; i<left.size(); i++)
    val(i) = CppAD::CondExpLt(left(i),right(i),left(i),right(i));
  return val;
}
template<class Type>
vector<Type> pmax(vector<Type> left, vector<Type> right) {
  vector<Type> val(left.size());
  for(int i=0; i<left.size(); i++)
    val(i) = CppAD::CondExpGt(left(i),right(i),left(i),right(i));
  return val;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(betas);
  DATA_MATRIX(X);
  DATA_MATRIX(XD);
  DATA_MATRIX(q_matrix);
  DATA_VECTOR(event);
  DATA_VECTOR(t);
  DATA_VECTOR(boundaryKnots);
  DATA_VECTOR(interiorKnots);
  ns<Type> s(boundaryKnots, interiorKnots, q_matrix, 1);
  vector<Type> eta = X * beta;
  vector<Type> etaD = XD * beta;
  vector<Type> logtstar = log(t) - eta;
  matrix<Type> Xs = s.basis(logtstar);
  matrix<Type> XDs = s.basis(logtstar,1);
  vector<Type> etas = Xs * betas;
  vector<Type> etaDs = XDs * betas;
  // fix bounds on etaDs
  vector<Type> eps(etaDs.size());
  eps.setZero();
  eps = eps+Type(1.0e-8);
  Type pen = dot(pmin<Type>(etaDs,eps), pmin<Type>(etaDs, eps));
  etaDs = pmax<Type>(etaDs, eps);
  // fix bounds on etaD
  pen += dot(pmin<Type>(Type(1)/t-etaD,eps), pmin<Type>(Type(1)/t-etaD,eps));
  etaD = Type(1)/t - pmax<Type>(Type(1)/t-etaD, eps);
  vector<Type> H = exp(etas);
  vector<Type> logh = etas + log(etaDs) + log(Type(1)/t - etaD);
  Type f = pen - (dot(logh,event) - H.sum());
  return f;
}

