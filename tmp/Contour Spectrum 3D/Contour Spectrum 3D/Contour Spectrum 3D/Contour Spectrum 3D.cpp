///////////////////////////////////////////////////////////////////
//This program computes entropy of tetrahedron mesh using the contour spectrum method
//5 datasets are available here as example
//
//Any questions, please contact Tiantian Xu at txu04@students.poly.edu
//
///////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "point.h"
#include <vector>
using std::vector;

//read datasets by id,1=blunt, 2=comb, 3=delta, 4=post, 5=spx
void file_reader(int i);
//compute spline function for cell #cid, results are stored in func[]
void local_func(int cid,double *func);
//sort array and store the arrangment of positions
void bubbleSort(double arr[],int n, int *pos);
//compute the area of triangle with vertex p1,p2 and p3
double triangle_area(double *p1,double *p2, double *p3);
//compute global entropy with scalar values divided into hisize cells
double glb_entropy(int hisize);
//compute the area of vertex #tgtid in cell #cid
double cmptl(int cid, int tgtid);
//basis B-spline function on knot t, result stored in f1,f2 and f3
void spline(double t[],double *f1,double *f2,double *f3);
//spline function (not reduced version), 6 knots, 3 control points, return 5 functions on five intervals respectively
void bspline(double t[],Point cp[],double *f1,double *f2,double *f3,double *f4,double *f5);
//efficient way to compute spline function; only compute related part of basis spline functions, return 3 functions on 3 intervals respectively
void efbspline(double t[],Point cp[],double *f1,double *f2,double *f3);
//integration of quadratic function f:y=a*x^2+b*x+c on [min,max]
double integr(double f[],double min, double max);
//linear interpolation using scalar value on segment p1p2
void linterp(double s1,double s2,double s,double p1[],double p2[],double *p);
//search for position of "value" in list[]
int binarysearch(double list[],double value, int low, int high,int len);
//compute contribution of a cell on scalar interval [min,max]
double contribution(double func[],double min, double max);

//define global variables
vector <int> cell_list(0);
vector <float> vertex_list(0);
int cn=0,vn=0;
//When efficient is 'true', use the efficient way to compute b-splines
bool efficient=true;

int _tmain(int argc, _TCHAR* argv[])
{ 
	   
	int i=0;
	double r;
	clock_t start,finish;
	double totaltime;

	for(i=1;i<6;i++)
	{
		cell_list.clear();
		vertex_list.clear();
		cn=0;
		vn=0;

	file_reader(i);

	printf("%d:\n",i);

	efficient=false;

	start=clock();

	r=glb_entropy(100);

	finish=clock();


	printf("Not Efficient: %f\n",r);

    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;

	r=r/(log(100.0)/log(2.0));
	printf("Divided by log(hisize): %f\n",r);
	

    printf("Running time: %f\n",totaltime);

	start=clock();

	r=glb_entropy(1000);

	finish=clock();


	printf("Not Efficient: %f\n",r);

    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;

    printf("Running time: %f\n",totaltime);
		r=r/(log(1000.0)/log(2.0));
	printf("Divided by log(hisize): %f\n",r);
	
	efficient=true;

	start=clock();

	r=glb_entropy(100);

	finish=clock();

	printf("Efficient: %f\n",r);

   
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
		r=r/(log(100.0)/log(2.0));
	printf("Divided by log(hisize): %f\n",r);
	
    printf("Running time: %f\n",totaltime);

	start=clock();

	r=glb_entropy(1000);

	finish=clock();

	printf("Efficient: %f\n",r);

   
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
		r=r/(log(1000.0)/log(2.0));
	printf("Divided by log(hisize): %f\n",r);
	
    printf("Running time: %f\n",totaltime);

	}
	
	return 0; 
}

void local_func(int cid,double *func)
{
	
   double areav2=cmptl(cid,1);

   int v1,v2,v3,v4;
	v1=cell_list[cid*4];
	v2=cell_list[cid*4+1];
	v3=cell_list[cid*4+2];
	v4=cell_list[cid*4+3];

	double list[4]={vertex_list[v1*4+3],vertex_list[v2*4+3],vertex_list[v3*4+3],vertex_list[v4*4+3]};
	int pos[4]={v1,v2,v3,v4};
	bubbleSort(list,4,pos);//pos stores the index of sorted elements, list[pos[0]]<list[pos[1]]<...

	//6 knots, 3 control points, degree 2
	double t[6]={0};

	t[0]=list[0];
	t[1]=(list[0]+list[1])/2;
	t[2]=list[1];
	t[3]=list[2];
	t[4]=(list[2]+list[3])/2;
	t[5]=list[3];

	double x=areav2/(list[1]-list[0]);
	
	double px=(list[1]+list[2])/2.0;
	double py=areav2*(list[2]-list[0])/(list[1]-list[0]);

	double value3=py*(list[3]-list[2])/(list[3]-list[1]);

    Point cp[3];
	cp[0].x=t[2];
	cp[0].y=areav2/2.0;
	cp[0].z=0;

	cp[1].x=px;
	cp[1].y=py;
	cp[1].z=0;

	cp[2].x=t[3];
	cp[2].y=value3/2.0;
	cp[2].z=0;

	if(!efficient){
	double f1[3],f2[3],f3[3],f4[3],f5[3];

	bspline(t,cp,f1,f2,f3,f4,f5);

	int i=0;

	for(i=0;i<6;i++)
	func[i]=t[i];
	
	for(i=6;i<9;i++)
    func[i]=f1[i-6];

	for(i=9;i<12;i++)
    func[i]=f2[i-9];

	for(i=12;i<15;i++)
    func[i]=f3[i-12];
				
	for(i=15;i<18;i++)
    func[i]=f4[i-15];
					
	for(i=18;i<21;i++)
    func[i]=f5[i-18];}
	else
	{
	double f1[3],f2[3],f3[3];

	efbspline(t,cp,f1,f2,f3);

	int i=0;

	for(i=0;i<6;i++)
	func[i]=t[i];
	
	for(i=6;i<9;i++)
    func[i]=f1[i-6];

	for(i=9;i<12;i++)
    func[i]=f2[i-9];

	for(i=12;i<15;i++)
    func[i]=f3[i-12];
	
	}

}

void bubbleSort(double arr[],int n, int *pos)
{ 
	double temp=0;
	int t=0;
	for(int i=0;i<n;i++)
		for(int j=0;j<n-1;j++)
		{
		  if(arr[j]>arr[j+1]) 
		  {
		    temp=arr[j];
			arr[j]=arr[j+1];
			arr[j+1]=temp;

			t=pos[j];
			pos[j]=pos[j+1];
			pos[j+1]=t;

		  }
		
		}
  

}

//calculate area of the triangle defined by p1,p2 and p3
double triangle_area(double *p1,double *p2, double *p3)
{
	double AB=(p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]);
	double BC=(p1[0]-p3[0])*(p1[0]-p3[0])+(p1[1]-p3[1])*(p1[1]-p3[1])+(p1[2]-p3[2])*(p1[2]-p3[2]);
	double CA=(p3[0]-p2[0])*(p3[0]-p2[0])+(p3[1]-p2[1])*(p3[1]-p2[1])+(p3[2]-p2[2])*(p3[2]-p2[2]);

	AB=sqrt(AB);
	BC=sqrt(BC);
	CA=sqrt(CA);

	double p=(AB+BC+CA)/2;

    double area=p*(p-AB)*(p-BC)*(p-CA);
	area=sqrt(area);

	return area;
}

void file_reader(int i)
{
	FILE *p;
	switch(i)
	{ 
	   case 1:
	p=fopen("blunt.scalar","r");
	break;
		case 2:
	p=fopen("comb.scalar","r");
	break;
		case 3:
	p=fopen("delta.scalar","r");
	break;
		case 4:
	p=fopen("post.scalar","r");
	break;
		case 5:
	p=fopen("spx.scalar","r");
	break;

	}

	if(p==NULL){ printf("NO SUCH FILE!"); //exit(0);
	}
	else{

	//cn=number of cells, vn=number of vertices
	fscanf(p,"%d%d",&vn,&cn);

	int i=0,j=0;
	float v;
		int c;

	for(i=0;i<vn;i++)
		for(j=0;j<4;j++)
		{	fscanf(p,"%f",&v);
			vertex_list.push_back(v);
			
	   }

	for(i=0;i<cn;i++)
		for(j=0;j<4;j++)
			{
				fscanf(p,"%d",&c);
				cell_list.push_back(c);
		}
	}

	fclose(p);
}

double glb_entropy(int hisize)
{
	
	int i=0,j=0;
	//extract scalar list


	double *scalar_list=new double[vn];
	for(i=0;i<vn;i++) 
		scalar_list[i]=vertex_list[i*4+3];

    double rangemin=scalar_list[0],rangemax=scalar_list[0];
	for(i=1;i<vn;i++)
	{
		if(scalar_list[i]>rangemax) rangemax=scalar_list[i];
		if(scalar_list[i]<rangemin) rangemin=scalar_list[i];
	}

    

	double list[1001]={0};
	for(i=0;i<hisize+1;i++)
	{	list[i]=rangemin+(rangemax-rangemin)*i/hisize;//hisize bins, bin i starts with scalar value list[i]
	}

	//process cell by cell: calculate the contribution of cell i to bin j
	double min=0,max=0;
	int posmin=0,posmax=0;
	
	double histogram[1001]={0};
	double func[21]={0};
	double tt=0;

    

	for(i=0;i<cn;i++)
	{
	   local_func(i,func);

	   //min and max defines range of cell i
       min=func[0];
	   max=func[5];

	   if(!efficient){

	   //cell i spans from bin[posmin] to bin[posmax], for every bin in this range, calculate corresponding contribution
	   posmin=binarysearch(list,min,0,hisize,hisize+1);
	   posmax=binarysearch(list,max,0,hisize,hisize+1);

	   for(j=posmin;j<=posmax;j++)
		  { 
			  tt=contribution(func,list[j],list[j+1]);
	       if(tt>=0||tt<=0)
			  histogram[j]+=tt;
	      }
	   }
	   else
	   {
	     posmin=binarysearch(list,min,0,hisize,hisize+1);

	     for(j=posmin;list[j]<max;j++)
		  { 
			  tt=contribution(func,list[j],list[j+1]);
	        if(tt>=0||tt<=0)
			  histogram[j]+=tt;
	      }
	   
	   }
	}

	double total=0;
	for(i=0;i<hisize;i++)
	{
	   total+=histogram[i];
	}

	double prb[1000]={0};
	for(i=0;i<hisize;i++)
		prb[i]=histogram[i]/total;

	double entropy=0;
		
	for(i=0;i<hisize;i++)
	{
		if(prb[i]>0)
	  entropy-=prb[i]*(log(prb[i])/log(2.0));
	}

	return entropy;
}

double cmptl(int cid, int tgtid)
{
	double area=0;

    //retrieve from cell list the vertices in current triangle
	int v1,v2,v3,v4;
	v1=cell_list[cid*4];
	v2=cell_list[cid*4+1];
	v3=cell_list[cid*4+2];
	v4=cell_list[cid*4+3];

	double list[4]={vertex_list[v1*4+3],vertex_list[v2*4+3],vertex_list[v3*4+3],vertex_list[v4*4+3]};
	int pos[4]={v1,v2,v3,v4};
	bubbleSort(list,4,pos);//pos stores the index of sorted elements, list[pos[0]]<list[pos[1]]<...


	
	double p1[3]={vertex_list[pos[0]*4],vertex_list[pos[0]*4+1],vertex_list[pos[0]*4+2]};
	double p2[3]={vertex_list[pos[1]*4],vertex_list[pos[1]*4+1],vertex_list[pos[1]*4+2]};
	double p3[3]={vertex_list[pos[2]*4],vertex_list[pos[2]*4+1],vertex_list[pos[2]*4+2]};
	double p4[3]={vertex_list[pos[3]*4],vertex_list[pos[3]*4+1],vertex_list[pos[3]*4+2]};

	double temp;
	double pi11[3],pi21[3];
	if(tgtid==1)
	{
	 linterp(list[0],list[2],list[1],p1,p3,pi11);
	 linterp(list[0],list[3],list[1],p1,p4,pi21);
	 area=triangle_area(pi11,pi21,p2);
	return area;
	
	}
	else if(tgtid==2)
	{
	 linterp(list[1],list[3],list[2],p2,p4,pi11);
	 linterp(list[0],list[3],list[2],p1,p4,pi21);
	 area=triangle_area(pi11,pi21,p3);
	 return area;
	}
    
	return area;
}

void spline(double t[],double *f1,double *f2,double *f3)
{

	double a=0,b=0,c=0;

	//part 1:t0<=t<t1
	double m=0,n=0;
	if(t[2]!=t[0])
	m=1/(t[2]-t[0]);
	if(t[1]!=t[0])
	n=1/(t[1]-t[0]);

	a=m*n;
	b=-2*t[0]*m*n;
	c=t[0]*t[0]*m*n;

	f1[0]=a;
	f1[1]=b;
	f1[2]=c;

	//part 2: t1<=t<t2
	double u=0,v=0,w=0;
	if(t[2]!=t[0])
	u=1/(t[2]-t[0]);
	if(t[2]!=t[1])
	v=1/(t[2]-t[1]);
	if(t[3]!=t[1])
	w=1/(t[3]-t[1]);

	a=-(u*v+w*v);
	b=u*v*(t[0]+t[2])+w*v*(t[1]+t[3]);
	c=-u*v*t[0]*t[2]-v*t[1]-w*v*t[1]*t[1];

	f2[0]=a;
	f2[1]=b;
	f2[2]=c;

	//part 3: t2<=t<t3
	m=0,n=0;
	if(t[3]!=t[1])
	m=1/(t[3]-t[1]);
	if(t[3]!=t[2])
	n=1/(t[3]-t[2]);

	a=m*n;
	b=-2*m*n*t[3];
	c=m*n*t[3]*t[3];

	f3[0]=a;
	f3[1]=b;
	f3[2]=c;

}

void bspline(double t[],Point cp[],double *f1,double *f2,double *f3,double *f4,double *f5)
{
	
	double knot1[4]={t[0],t[1],t[2],t[3]};
	double knot2[4]={t[1],t[2],t[3],t[4]};
	double knot3[4]={t[2],t[3],t[4],t[5]};

	double f11[3],f12[3],f13[3];
	double f21[3],f22[3],f23[3];
	double f31[3],f32[3],f33[3];

	spline(knot1,f11,f12,f13);
	spline(knot2,f21,f22,f23);
	spline(knot3,f31,f32,f33);


	//(t0,t1):  y=f11*cp[0].y
	f1[0]=f11[0]*cp[0].y;
	f1[1]=f11[1]*cp[0].y;
	f1[2]=f11[2]*cp[0].y;	

	//(t1,t2): y=f12*cp[0].y+f21*cp[1].y
	f2[0]=f12[0]*cp[0].y+f21[0]*cp[1].y;
	f2[1]=f12[1]*cp[0].y+f21[1]*cp[1].y;
	f2[2]=f12[2]*cp[0].y+f21[2]*cp[1].y;

	//(t2,t3): y=f13*cp[0].y+f22*cp[1].y+f31*cp[2].y
	f3[0]=f13[0]*cp[0].y+f22[0]*cp[1].y+f31[0]*cp[2].y;
	f3[1]=f13[1]*cp[0].y+f22[1]*cp[1].y+f31[1]*cp[2].y;
	f3[2]=f13[2]*cp[0].y+f22[2]*cp[1].y+f31[2]*cp[2].y;

	//(t3,t4): y=f23*cp[1].y+f32*cp[2].y
	f4[0]=f23[0]*cp[1].y+f32[0]*cp[2].y;
	f4[1]=f23[1]*cp[1].y+f32[1]*cp[2].y;
	f4[2]=f23[2]*cp[1].y+f32[2]*cp[2].y;

	//(t4,t5): y=f33*cp[2].y
	f5[0]=f33[0]*cp[2].y;
	f5[1]=f33[1]*cp[2].y;
	f5[2]=f33[2]*cp[2].y;
       

}

void efbspline(double t[],Point cp[],double *f1,double *f2,double *f3)
{
	double b11[3],b13[3],b22[3],b31[3],b33[3];
	double a=0,b=0,c=0;
	double m=0,n=0;

	//for basis b-spline function on t0-t3,only compute first and third part
	//part 1:t0<=t<t1
	
	if(t[2]!=t[0])
	m=1/(t[2]-t[0]);
	if(t[1]!=t[0])
	n=1/(t[1]-t[0]);

	a=m*n;
	b=-2*t[0]*m*n;
	c=t[0]*t[0]*m*n;

	b11[0]=a;
	b11[1]=b;
	b11[2]=c;


	//part 3: t2<=t<t3
	m=0,n=0;
	if(t[3]!=t[1])
	m=1/(t[3]-t[1]);
	if(t[3]!=t[2])
	n=1/(t[3]-t[2]);

	a=m*n;
	b=-2*m*n*t[3];
	c=m*n*t[3]*t[3];

	b13[0]=a;
	b13[1]=b;
	b13[2]=c;

	//for basis b-spline function on t1-t4, compute the second part
	//t2<=t<=t3
	double u=0,v=0,w=0;
	if(t[3]!=t[1])
	u=1/(t[3]-t[1]);
	if(t[3]!=t[2])
	v=1/(t[3]-t[2]);
	if(t[4]!=t[2])
	w=1/(t[4]-t[2]);

	a=-(u*v+w*v);
	b=u*v*(t[1]+t[3])+w*v*(t[2]+t[4]);
	c=-u*v*t[1]*t[3]-v*t[2]-w*v*t[2]*t[2];

	b22[0]=a;
	b22[1]=b;
	b22[2]=c;

	//for basis b-spline function on t2-t5, compute the first and third part
	//part 1:t2<=t<t3
	
	if(t[4]!=t[2])
	m=1/(t[4]-t[2]);
	if(t[3]!=t[2])
	n=1/(t[3]-t[2]);

	a=m*n;
	b=-2*t[2]*m*n;
	c=t[2]*t[2]*m*n;

	b31[0]=a;
	b31[1]=b;
	b31[2]=c;


	//part 3: t4<=t<t5
	m=0,n=0;
	if(t[5]!=t[3])
	m=1/(t[5]-t[3]);
	if(t[5]!=t[4])
	n=1/(t[5]-t[4]);

	a=m*n;
	b=-2*m*n*t[5];
	c=m*n*t[5]*t[5];

	b33[0]=a;
	b33[1]=b;
	b33[2]=c;

	//compute b-spline function
	//(t0,t2): y=b11*cp[0].y
	f1[0]=b11[0]*cp[0].y;
	f1[1]=b11[1]*cp[0].y;
	f1[2]=b11[2]*cp[0].y;	
	
	//(t2,t3): y=b13*cp[0].y+b22*cp[1].y+b31*cp[2].y
	f2[0]=b13[0]*cp[0].y+b22[0]*cp[1].y+b31[0]*cp[2].y;
	f2[1]=b13[1]*cp[0].y+b22[1]*cp[1].y+b31[1]*cp[2].y;
	f2[2]=b13[2]*cp[0].y+b22[2]*cp[1].y+b31[2]*cp[2].y;

	//(t3,t5): y=b33*cp[2].y
	f3[0]=b33[0]*cp[2].y;
	f3[1]=b33[1]*cp[2].y;
	f3[2]=b33[2]*cp[2].y;



}

int binarysearch(double list[],double value,int low, int high,int len)
{
	  if(value==list[high]) return high;
	
      int mid = low + (high - low) / 2;
            if (list[0] > value)
            {
                return -1;
            }
            else if (list[len-1] <= value)
            {
			
                return  - 1;
            }
               
            else if (list[mid] > value)
            {
                if (list[mid - 1] < value)
                {
                    return mid - 1;
                }
                else
                    return binarysearch(list, value, low, mid - 1,len);
            }

            else if (list[mid] < value && !(list[mid + 1] > value))
            {
                return binarysearch(list, value, mid + 1, high,len);
            }

            else return mid;

}

double contribution(double func[],double min, double max)
{
	if(!efficient){
	int i=0,p=0;
	double v[6],f1[3],f2[3],f3[3],f4[3],f5[3];
	for(i=0;i<6;i++) v[i]=func[i];
	for(i=0;i<3;i++) f1[i]=func[i+6];
	for(i=0;i<3;i++) f2[i]=func[i+9];
	for(i=0;i<3;i++) f3[i]=func[i+12];
	for(i=0;i<3;i++) f4[i]=func[i+15];
	for(i=0;i<3;i++) f5[i]=func[i+18];

	double result=0;
	
    if(min<=v[0])
	{ 
	  if(max<=v[0])
		  return 0;
	  else if(max<=v[1])
		  result+=integr(f1,v[0],max);
	  else if(max<=v[2])
	  {
	     result+=integr(f1,v[0],v[1]);
		 result+=integr(f2,v[1],max);
	  }else if(max<=v[3])
	  {
	     result+=integr(f1,v[0],v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],max);
	  
	  }else if(max<=v[4])
	  {
	     result+=integr(f1,v[0],v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],max);
	  }else if(max<=v[5])
	  {
	     result+=integr(f1,v[0],v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],v[4]);
		 result+=integr(f5,v[4],max);
	  }else if(max>v[5])
	  {
	     result+=integr(f1,v[0],v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],v[4]);
		 result+=integr(f5,v[4],v[5]);
	  
	  }else return 0;
	  return result;
	}
	else if(min<=v[1])
	{
	    if(max<=v[1])
		  result+=integr(f1,min,max);
	  else if(max<=v[2])
	  {
	     result+=integr(f1,min,v[1]);
		 result+=integr(f2,v[1],max);
	  }else if(max<=v[3])
	  {
	     result+=integr(f1,min,v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],max);
	  
	  }else if(max<=v[4])
	  {
	     result+=integr(f1,min,v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],max);
	  }else if(max<=v[5])
	  {
	     result+=integr(f1,min,v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],v[4]);
		 result+=integr(f5,v[4],max);
	  }else if(max>v[5])
	  {
	     result+=integr(f1,min,v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],v[4]);
		 result+=integr(f5,v[4],v[5]);
	  
	  }else return 0;
	  
	  return result;
	
	}
	else if(min<=v[2])
	{
	   if(max<=v[2])
	  {
		 result+=integr(f2,min,max);
	  }else if(max<=v[3])
	  {
		 result+=integr(f2,min,v[2]);
		 result+=integr(f3,v[2],max);
	  
	  }else if(max<=v[4])
	  {
		 result+=integr(f2,min,v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],max);
	  }else if(max<=v[5])
	  {
		 result+=integr(f2,min,v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],v[4]);
		 result+=integr(f5,v[4],max);
	  }else if(max>v[5])
	  {
		 result+=integr(f2,min,v[2]);
		 result+=integr(f3,v[2],v[3]);
		 result+=integr(f4,v[3],v[4]);
		 result+=integr(f5,v[4],v[5]);
	  
	  }else return 0;
	  return result;
	}
	else if(min<=v[3])
	{
	   if(max<=v[3])
	  {
		 result+=integr(f3,min,max);
	  
	  }else if(max<=v[4])
	  {
		 result+=integr(f3,min,v[3]);
		 result+=integr(f4,v[3],max);
	  }else if(max<=v[5])
	  {
		 result+=integr(f3,min,v[3]);
		 result+=integr(f4,v[3],v[4]);
		 result+=integr(f5,v[4],max);
	  }else if(max>v[5])
	  {
		 result+=integr(f3,min,v[3]);
		 result+=integr(f4,v[3],v[4]);
		 result+=integr(f5,v[4],v[5]);
	  
	  }else return 0;

	  return result;
	}
	else if(min<=v[4])
	{
		if(max<=v[4])
	  {
		 result+=integr(f4,min,max);
	  }else if(max<=v[5])
	  {
		 result+=integr(f4,min,v[4]);
		 result+=integr(f5,v[4],max);
	  }else if(max>v[5])
	  {
		 result+=integr(f4,min,v[4]);
		 result+=integr(f5,v[4],v[5]);
	  
	  }else return 0;

	  return result;
	}
	else if(min<v[5])
	{
	  if(max<=v[5]) result+=integr(f5,min,max);
	  else if(max>v[5]) result+=integr(f5,min,v[5]);
	  else return 0;

	  return result;
	}
	else return 0;

    }
else{
	int i=0,p=0;
	double v[4],f1[3],f2[3],f3[3];
	double result=0;
	v[0]=func[0];
	v[1]=func[2];
	v[2]=func[3];
	v[3]=func[5];

	for(i=0;i<3;i++) f1[i]=func[i+6];
	for(i=0;i<3;i++) f2[i]=func[i+9];
	for(i=0;i<3;i++) f3[i]=func[i+12];

	if(min<=v[0])
	{
	  if(max<=v[0]) return 0;
	  else if(max<=v[1])
		  result+=integr(f1,v[0],max);
	  else if(max<=v[2])
	  {
	     result+=integr(f1,v[0],v[1]);
		 result+=integr(f2,v[1],max);
	  }else if(max<=v[3])
	  {
	     result+=integr(f1,v[0],v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],max);
	  
	  }else if(max>v[3])
	  {
	     result+=integr(f1,v[0],v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],v[3]);
	  
	  }else return 0;
	  return result;
	}
	else if(min<=v[1])
	{
	    if(max<=v[1])
		  result+=integr(f1,min,max);
	  else if(max<=v[2])
	  {
	     result+=integr(f1,min,v[1]);
		 result+=integr(f2,v[1],max);
	  }else if(max<=v[3])
	  {
	     result+=integr(f1,min,v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],max);
	  
	  }else if(max>v[3])
	  {
	     result+=integr(f1,min,v[1]);
		 result+=integr(f2,v[1],v[2]);
		 result+=integr(f3,v[2],v[3]);
	  
	  }else return 0;
	
	return result;
	}
	else if(min<=v[2])
	{
	   if(max<=v[2])
	  {
		 result+=integr(f2,min,max);
	  }else if(max<=v[3])
	  {
		 result+=integr(f2,min,v[2]);
		 result+=integr(f3,v[2],max);
	  
	  }else if(max>v[3])
	  {
		 result+=integr(f2,min,v[2]);
		 result+=integr(f3,v[2],v[3]);
	  
	  }else return 0;
	 return result;
	}
	else if(min<=v[3])
	{
	   if(max<=v[3])
	  {
		 result+=integr(f3,min,max);
	  
	  }else  if(max>v[3])
	  {
		 result+=integr(f3,min,v[3]);
	  
	  }else return 0;
	  return result;
	}else return 0;
	
	}
	 
}

double integr(double f[],double min, double max)
{
	if(min==max||min>max) return 0;

		double result=0;

	double p3=0,p2=0,p1=0;

	p3=f[0]/3; 
	p2=f[1]/2;
	p1=f[2];
	
	result=p3*(pow(max,3)-pow(min,3))+p2*(pow(max,2)-pow(min,2))+p1*(max-min);

	return result;


}

void linterp(double s1,double s2,double s,double p1[],double p2[],double *p)
{
     double ratio=(s-s1)/(s2-s1);

	 for(int i=0;i<3;i++)
		 p[i]=p1[i]+ratio*(p2[i]-p1[i]);


}