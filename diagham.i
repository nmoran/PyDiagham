%module diagham
%include cpointer.i 

%{
#define SWIG_FILE_WITH_INIT
#include "RealVector.h"
#include "RealMatrix.h"
#include "BosonOnSphereTwoLandauLevels.h"
#include "utils.h"
%}


ostream& GetCout();

class RealVector
{
 public:

  RealVector();
  RealVector(int size, bool zeroFlag = false);  
  RealVector(double* array, int size);
  RealVector(const RealVector& vector, bool duplicateFlag = false);  
  RealVector(const Vector& vector);
  ~RealVector ();
  void Resize (int dimension);
  void Resize (long dimension);
  void ResizeAndClean (int dimension);
  RealVector& Normalize();
  RealVector Extract(int firstCoordinate, int lastCoordinate, int step = 1);
  bool WriteVector (const char* fileName);
  bool WriteAsciiVector (const char* fileName);
  bool ReadVector (const char* fileName);
  bool ReadVector (const char* fileName, long minIndex, long maxIndex);
  long ReadVectorDimension (const char* fileName);
  bool ReadVectorTest (const char* fileName);
};


class RealMatrix 
{
 public:
	RealMatrix();
	RealMatrix(int nbrRow, int nbrColumn, bool zero = false);
	RealMatrix(RealVector* columns, int nbrColumn);
	RealMatrix(const RealMatrix& M);
	RealMatrix(Matrix& M);
	~RealMatrix();
	void GetMatrixElement(int i, int j, double& x) const;
	void SetMatrixElement(int i, int j, double x);
	void SetMatrixElement(int i, int j, const Complex& x);
	void AddToMatrixElement(int i, int j, double x);
	void AddToMatrixElement(int i, int j, const Complex& x);
	void Resize (int nbrRow, int nbrColumn);
	void ResizeAndClean (int nbrRow, int nbrColumn);
	RealMatrix& NormalizeColumns ();
	RealMatrix& Transpose ();
	double Determinant ();
	double Permanent();
	double* SingularValueDecomposition(RealMatrix& uMatrix, RealMatrix& vMatrix);
	double* SingularValueDecomposition();
};


class BosonOnSphereTwoLandauLevels{
  public:
	BosonOnSphereTwoLandauLevels (int nbrBosons, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);
	~BosonOnSphereTwoLandauLevels ();
	int GetParticleStatistic();
	ostream& PrintState ( ostream& Str, int state);
	ostream& PrintStateMonomial ( ostream& Str, int state);
	int GetHilbertSpaceDimension();
	RealMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, bool removeBinomialCoefficient);
};

class FermionOnSphereTwoLandauLevels{
  public:
    FermionOnSphereTwoLandauLevels (int nbrBosons, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);
	~FermionOnSphereTwoLandauLevels ();
	int GetParticleStatistic();
	ostream& PrintState ( ostream& Str, int state);
	ostream& PrintStateMonomial ( ostream& Str, int state);
	int GetHilbertSpaceDimension();
	double AddAd (int index, int m);
	double AduAu (int index, int m);
	int AduAu (int index, int m, int n, double& coefficient);
	int AddAd (int index, int m, int n, double& coefficient);
	int AduAd (int index, int m, int n, double& coefficient);
	int AddAu (int index, int m, int n, double& coefficient);
	double AuAu (int index, int n1, int n2);
	double AdAd (int index, int n1, int n2);
	double AuAd (int index, int n1, int n2);
	int AduAdu (int m1, int m2, double& coefficient);
	int AddAdd (int m1, int m2, double& coefficient);
	int AduAdd (int m1, int m2, double& coefficient);
};


/*class ParticleOnSphereTwoLandauLevelDeltaHamiltonian{
  public:
	ParticleOnSphereTwoLandauLevelDeltaHamiltonian (ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
							double** pseudoPotential, double *cyclotronEnergy,
							AbstractArchitecture* architecture,long memory, bool onDiskCacheFlag, char* precalculationFileName);
	~ParticleOnSphereTwoLandauLevelDeltaHamiltonian ();
	void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination);
	int GetParticleStatistic();
	#        void PrintStateSimple (int state);
	#        void PrintStateMonomialSimple ( int state);
	ostream& PrintState ( ostream& Str, int state);
	ostream& PrintStateMonomial ( ostream& Str, int state);
	int GetHilbertSpaceDimension();
};*/
