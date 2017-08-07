#ifndef _MCpoint_h_
#define _MCpoint_h_

#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <mex.h>
#include <szMexUtility.h>

const int MCpointMaxCandidates=10;

class MCpoint
{
public:
   MCpoint(float y=0, float x=0, float vy=0, float vx=0, int nump=1, int numc=1) 
   {
      Uy=vy;
      Ux=vx;
      X=x;
      Y=y;
      NumPoints=nump;
      NumCandidates=numc;
      Vy=vector<float>(numc,0);
      Vx=vector<float>(numc,0);
      Prob=vector<float>(numc+1,1.0/(numc+1));
	  Index = vector<int>(numc);
      int i;
      int numlink=MCpointMaxCandidates*MCpointMaxCandidates*nump;
      Link=vector<float>(numlink,.5);
      for(i=0; i<numlink; ++i)
         Link[i]=0.5;
   }

   /*MCpoint(const MCpoint& cp) 
   {
      Ux=cp.GetVelocityY();
      Uy=cp.GetVelocityX();
      X=cp.GetX();
      Y=cp.GetY();
      
      NumPoints=cp.GetNumPoints();
      NumCandidates=cp.GetNumCandidates();
      assert(NumCandidates<=MCpointMaxCandidates);
      Vy=vector<float>(NumCandidates,0);
      Vx=vector<float>(NumCandidates,0);
      Prob=vector<float>(NumCandidates+1,1.0/(NumCandidates+1));
      Index=vector<int>(NumCandidates);
      int i;
      for(i=0; i<NumCandidates; ++i) 
	  {
         Vy[i]=cp.GetEstimateVelocityY(i);
         Vx[i]=cp.GetEstimateVelocityX(i);
         Index[i]=cp.GetIndex(i);
      }
	  for (i = 0; i <= NumCandidates; ++i)
	  {
		  Prob[i] = cp.GetProb(i);
	  }
      int numlink=MCpointMaxCandidates*MCpointMaxCandidates*NumPoints;
      Link=vector<float>(numlink);
	  for (i = 0; i < numlink; ++i)
	  {
		  SetLink(i, cp.GetLink(i));
	  }
   }*/
   
   /*const MCpoint& operator=(const MCpoint& cp) 
   {
      Ux=cp.GetVelocityY();
      Uy=cp.GetVelocityX();
      X=cp.GetX();
      Y=cp.GetY();
      
      NumPoints=cp.GetNumPoints();
      NumCandidates=cp.GetNumCandidates();
      assert(NumCandidates<=MCpointMaxCandidates);
      Vy=vector<float>(NumCandidates,0);
      Vx=vector<float>(NumCandidates,0);
      Prob=vector<float>(NumCandidates+1,1.0/(NumCandidates+1));
	  Index = vector<int>(NumCandidates);
      int i;
      for(i=0; i<NumCandidates; ++i) 
	  {
         Vy[i]=cp.GetEstimateVelocityY(i);
         Vx[i]=cp.GetEstimateVelocityX(i);
         Index[i]=cp.GetIndex(i);
      }
      for(i=0; i<=NumCandidates; ++i)
         Prob[i]=cp.GetProb(i);
      
      int numlink=MCpointMaxCandidates*MCpointMaxCandidates*NumPoints;
      Link=vector<float>(numlink);
      for(i=0; i<numlink; ++i)
         SetLink(i,cp.GetLink(i));
      
      return *this;
   }*/
   
   /*~MCpoint() {
   }*/
   
   
   inline float GetVelocityY() const {return Uy;}
   inline float GetVelocityX() const {return Ux;}
   inline float GetX() const {return X;}
   inline float GetY() const {return Y;}
   inline int GetNumPoints() const {return NumPoints;}
   inline int GetNumCandidates() const {return NumCandidates;}
   inline float GetLink(int j, int m, int n) const 
   {
      return Link[j*MCpointMaxCandidates*MCpointMaxCandidates+m*MCpointMaxCandidates+n];
   }
   inline float GetLink(int n) const {return Link[n];}
   inline float GetProb(int i) const {return Prob[i];}
   inline int GetIndex(int i) const {return Index[i];}
   inline float GetEstimateVelocityY(int i) const {return Vy[i];}
   inline float GetEstimateVelocityX(int i) const {return Vx[i];}
   pair<float, int> GetDomiantCandidate() {
	   pair<float, int> result(0.0f, -1);
	   for (int i = 0; i < Prob.size(); ++i)
	   {
		   if (result.first < Prob[i])
		   {
			   result.first = Prob[i];
			   result.second = i < NumCandidates ? Index[i] : -1;
		   }
	   }
	   return result;
   }
   
   inline void SetEstimateVelocityY(int n, float v) 
   {
      Vy[n]=v;
   }
   inline void SetEstimateVelocityX(int n, float v) 
   {
      Vx[n]=v;
   }
   inline void SetProb(int n, float p) 
   {
      Prob[n]=p;
   }
   inline void SetIndex(int m, int j) 
   {
      Index[m]=j;
   }
   inline void SetVelocityY(float v) 
   {
      Uy=v;
   }
   inline void SetVelocityX(float v)
   {
      Ux=v;
   }
   inline void SetX(float x) {
      X=x;
   }
   inline void SetY(float y) {
      Y=y;
   }

   inline void SetNumPoints(int n) 
   {
      if(NumPoints==n)
         return;
      NumPoints=n;
      int numlink=MCpointMaxCandidates*MCpointMaxCandidates*NumPoints;
      Link=vector<float>(numlink,.5);
   }
   
   inline void SetNumCandidates(int n)
   {
      if(NumCandidates==n)
         return;
      NumCandidates=n;
      Prob=vector<float>(NumCandidates+1);
      Vy=vector<float>(NumCandidates);
      Vx=vector<float>(NumCandidates);
	  Index = vector<int>(NumCandidates);
      int numlink=MCpointMaxCandidates*MCpointMaxCandidates*NumPoints;
      Link=vector<float>(numlink,.5);
   }
   
   inline void SetLink(int j, int m, int n, float link)
   {
      Link[j*MCpointMaxCandidates*MCpointMaxCandidates+m*MCpointMaxCandidates+n]=link;
   }
   
   inline void SetLink(int n, float link)
   {
      Link[n]=link;
   }
   
   inline void Move(float rate, float var)
   {
      Y+=rate*(Uy+var*(rndm(0)-.5));
      X+=rate*(Ux+var*(rndm(0)-.5));
   }

   void print() 
   {
	   printf("X=%f, Y=%f\n", GetX(), GetY());
	   printf("vX=%f, vY=%f\n", GetVelocityX(), GetVelocityY());
	   int i, j, k;
	   /*printf("Link Weight:\n");
	   for (i = 0; i < GetNumPoints(); ++i)
	   {
		   for (j = 0; j < GetNumCandidates(); ++j)
		   {

			   for (k = 0; k < MCpointMaxCandidates; ++k)
			   {
				   printf("%f, ", GetLink(i, j, k));
			   }
		   }
		   //printf(": %f\n", GetLink(i, 0, 0));
	   }*/
	   /*printf("Link: ");
	   for (int i = 0; i < Link.size(); ++i)
	   {
		   if (Abs(Link[i] - 0.5) > 0.01)
		   {
			   printf("%f, ", Link[i]);
		   }
	   }
	   printf("\n");*/
	   printf("Velocity Estimates:\n");
	   for (i = 0; i<GetNumCandidates(); ++i) 
	   {
		   printf("%f, %f, %f, %d\n", GetEstimateVelocityX(i), GetEstimateVelocityY(i), GetProb(i), GetIndex(i));
	   }
	   printf("Unmatched: %f\n", GetProb(i));
   }
   
 protected:
    vector<float> Vx; // velocity estimate
    vector<float> Vy;
    vector<float> Prob;
    vector<float> Link;
	vector<int> Index;
    float X;
    float Y;
    float Ux;  // true velocity
    float Uy;
    int NumPoints;
    int NumCandidates;
};

typedef vector<MCpoint> vMCpoint;

ostream& operator<<(ostream& s, const MCpoint& v) {
   s << "X=" << v.GetX() << ", Y=" << v.GetY();
   s << ", Vx=" << v.GetVelocityX() << ", Vy=" << v.GetVelocityY() << endl;
   int i,j,k;
   /*s << "Link Weight:" << endl;
   for(i=0; i<v.GetNumPoints(); ++i) {
      for(j=0; j<v.GetNumCandidates(); ++j)
         for(k=0; k<MCpointMaxCandidates; ++k)
            s << v.GetLink(i,j,k) << ", ";
         //s << v.GetLink(i,0,0) << ", ";
         s << endl;
   }*/
   s << "Velocity Estimates:" << endl;
   for(i=0; i<v.GetNumCandidates(); ++i) {
      s << v.GetEstimateVelocityX(i) << ", ";
      s << v.GetEstimateVelocityY(i) << ", ";
      s << v.GetProb(i) << ", ";
      s << v.GetIndex(i) << "\n";
   }
   s << endl;
   return s;
}
#endif /* MCpoint_h */
