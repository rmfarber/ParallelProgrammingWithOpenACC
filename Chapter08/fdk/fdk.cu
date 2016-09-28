/*  
******************************************************
This program is to reconstruct for 3-D cone beam projection, apply on 3-D shep-Logan head phaton
There are three steps to the weighted filtered backprojection algorithm: 
1) convert projection to projection_prime (weighted)
2) filtering part
3) backprojection part

reference book: "Principles of Computerized Tomographic Imaging"
                 Avinash C. Kak  & Malcolm Slaney    Page 100-107
implement in Frequency Domain

Permission to use, copy, distribute and modify this software for any 
purpose with or without fee is hereby granted. This software is        
provided "as is" without express or implied warranty. 

Send comments or suggestions for this OpenACC version to
            rxu6@uh.edu, schandra@udel.edu

Authors: Rengan Xu, Sunita Chandrasekaran

May 26th, 2016
******************************************************
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<string.h>
#include<cuda_runtime.h>

#define PI 3.141592653589793
#define SO 62.0  // distance from source to rotation center
#define SD 83.0 // distance from detector to source center. originally 83
#define OD (SO/SD)  // for translate real detector to the dector at origin
#define PROJECTION_Y  200    // number of projections in y axis of detector
#define PROJECTION_Z  200 // number of projections in z axis of detector
#define frame_size (PROJECTION_Z*PROJECTION_Y)
#define Z_CENTER 0  // 115.5
#define Y_CENTER 100 // 515.5
#define sample_interval_y (0.0194*OD)   // sample interval in image detector at origin (not true interval at real detector)
#define sample_interval_z (0.0194*OD)
#define ZP 256  // zero padding
#define RECONSIZE 200
#define RECONSIZE_Z 200
#define REC_XY_CENTER ((RECONSIZE-1.0)/2.0)
#define CONVOLVESIZE ZP/2   // CONVOLVESIZE = ZP/2
#define num_belta 300
#define belta_step (360.0/300.0) //1.0
#define recon_step (sample_interval_y) // voxel size
#define recon_step_z   recon_step
#define water 0.2006
#define ignore1 0
#define open 9013
#define zstart 200  // go up
#define THREADS 128

void four1(float data[], unsigned long nn, int isign);
void realft(float data[], unsigned long n, int isign);
void cosft1(float y[], int n);
//version 2.0
/************************************ Kernel for the 3rd step back projeciton *******************************/
__global__ void back_projection(float *fp_d, short int* CT_numbers_d, int view_start,
					 int view_end, int X_SIZE, int Y_SIZE, int Z_SIZE)
{
	float x,y,z,t,s,p_prime, ksi_prime,SO_s,factor,belta_rad,cos_belta,sin_belta;
	float m_f,m_z,n_f,n_y,temp;
	int m_less,n_less, i, j ,k, l, idx;
	short int CT_number;
	float rec;
    int size2;
    
    idx = threadIdx.x + blockIdx.x*blockDim.x;
    size2 = RECONSIZE*RECONSIZE;
//    size1 = RECONSIZE_Z*size2;

//for(idx=idx; idx<size1; idx+=blockDim.x*gridDim.x)
{
    i = (idx/size2)%RECONSIZE_Z;
    j = (idx/RECONSIZE)%RECONSIZE;
    k = idx%RECONSIZE;
	
    z=(Z_CENTER-zstart+i)*recon_step_z;
	y=(j-REC_XY_CENTER)*recon_step;
	x=(k-REC_XY_CENTER)*recon_step;
	rec = 0;
	for(l=view_start;l<view_end;l++)  
	{     
		belta_rad=(num_belta-l)*belta_step*PI/180;
		cos_belta=cos(belta_rad); sin_belta=sin(belta_rad);
		t=x*cos_belta+y*sin_belta;
		s=y*cos_belta-x*sin_belta;
		SO_s=SO/(SO-s);
		p_prime=SO_s*t;
		ksi_prime=SO_s*z;
		factor=SO_s*SO_s;
		/*    bilinear interpolation  */
		m_f=-ksi_prime/sample_interval_z+Z_CENTER;
		m_less=(int)floor(m_f);
		m_z=(m_f-m_less);
		n_f=p_prime/sample_interval_y+Y_CENTER;
		n_less=(int)floor(n_f);
		n_y=(n_f-n_less);
					
        if (m_less>=199) m_less = PROJECTION_Z-2;
        if (n_less>=255) n_less = ZP-2;
        if (m_less<=0) m_less = 0;
        if (n_less<=0) n_less = 0;

		temp=(1-m_z)*(1-n_y)*fp_d[l*Y_SIZE*Z_SIZE+m_less*Z_SIZE+n_less]+m_z*(1-n_y)*fp_d[l*Y_SIZE*Z_SIZE+(m_less+1)*Z_SIZE+n_less]+
			(1-m_z)*n_y*fp_d[l*Y_SIZE*Z_SIZE+m_less*Z_SIZE+n_less+1]+m_z*n_y*fp_d[l*Y_SIZE*Z_SIZE+(m_less+1)*Z_SIZE+n_less+1];

		rec+=factor*temp;
	}   // end of belta--viewend
	temp=rec*4*PI/num_belta;
	if(temp<0) 
		temp=0.0;
	CT_number=(short int)((temp-water)/water*1000);
	CT_numbers_d[i*RECONSIZE*RECONSIZE + j*RECONSIZE + k] = CT_number;
  }
}

/*************************************               MAIN    *************************************************/

int main(int argc, char** argv)  
{
    if(argc < 3)
    {
        printf("usage: %s <input> <output>\n", argv[0]);
        exit(0);
    }
	
    FILE *ptr_proj,*ptr_ct;

	int i,j,k,l,n,skip;
	static unsigned short proj;
	static float projection[ZP];
	static float weight[PROJECTION_Z][PROJECTION_Y];

    // y_prime & z_prime are coordinate in detector; y[] & z[] are coordinate in detector of rotation center
	float y_prime,z_prime,filter[CONVOLVESIZE+1];  

	int num_view,view_start,view_end;

	float *fp_h;
	float *fp_d;
	short int *CT_numbers_h;
	short int *CT_numbers_d;
	float ***filteredprojection;
    
    struct timeval tim;
    double begin, end;
    

	int X_SIZE = 300;
	int Y_SIZE = PROJECTION_Z;
	int Z_SIZE = ZP;
	filteredprojection = (float ***)malloc(sizeof(float **) * X_SIZE);
	 
	for (i = 0 ;  i < X_SIZE; i++) {
	   filteredprojection[i] = (float **)malloc(sizeof(float *) * Y_SIZE);
	 
	   for (j = 0; j < Y_SIZE; j++)
		  filteredprojection[i][j] = (float *)malloc(sizeof(float) * Z_SIZE);
	}

	// ramp filter design
	for(n=0;n<CONVOLVESIZE+1;n++)  
	{
		if(n==0) 
			filter[n]=1/(8*(sample_interval_y)*(sample_interval_y));
		else  
			if((n%2)==0) 
				filter[n]=0;
			else 
				filter[n]=-1/(2*n*n*PI*PI*(sample_interval_y)*(sample_interval_y));
	}
	cosft1(filter-1,CONVOLVESIZE);  // FFT

	for(i=0;i<PROJECTION_Z;i++)  
	{   // weitht factor is independent of belta rotation angle
		z_prime=-(i-Z_CENTER)*sample_interval_z;  // z center 414
		for(j=0;j<PROJECTION_Y;j++)  
		{
			y_prime=(j-Y_CENTER)*sample_interval_y;   // y center 504
			weight[i][j]=SO/sqrt(SO*SO+ y_prime *y_prime+z_prime*z_prime);
		}
	}
	
	num_view=num_belta;//process_size;
	view_start=0; //my_rank * num_view;
	view_end=view_start+num_view;

	if((ptr_proj=fopen(argv[1],"rb"))==NULL)  
	{ // If file open is not succesful, print could not open the file and quit
		fprintf(stderr,"Sorry could not open the file %s.\n", argv[1]);
		exit(1);
	}

	if((ptr_ct=fopen(argv[2],"wb"))==NULL ) 
	{ // If file open is not succesful, print could not open the file and quit
		fprintf(stderr,"Sorry could not open the file %s.\n", argv[2]);
		exit(1);
	}
    
    gettimeofday(&tim, NULL);
    begin = tim.tv_sec + (tim.tv_usec/1000000.0);

	// filtering projection data of each theta angle  in frequency domain
	for(l=view_start;l<view_end;l++)  
	{  // start of l ----------------
		skip=num_belta-l+ignore1;
		fseek(ptr_proj,(skip*frame_size)*sizeof(unsigned short),0);  // skip seq header , start from current frame l

		//step 1: convert projection to projection_prime  page 106 function (175)
		for(i=0;i<PROJECTION_Z;i++)  
		{
			for(j=0;j<ZP;j++)  
			{
				if( j<PROJECTION_Y)  
				{
					fread(&proj,sizeof(unsigned short),1,ptr_proj);
					//  if(j<50 || j>1900) proj=open; //correction for the collimators
					if( proj==0 )  
						proj=open;
					if(proj>=open) 
						projection[j]=0;
					else 
						projection[j]=-(log(proj*1.0/open))*weight[i][j];
				}
				else 
					projection[j]=0;    // zero padding
			}  //end of j-PROJECTION_Y
			realft(projection-1,ZP,1);       //  FFT

		   //step 2: filter projection
			for(j=0;j<ZP;j++)  
			{  // filter process in freq domain == filter * projection
				if((j%2)==0)  
					filteredprojection[l][i][j]=projection[j]*filter[j/2]*(sample_interval_y)*2/ZP;
				else  
					filteredprojection[l][i][j]=projection[j]*filter[(j-1)/2]*(sample_interval_y)*2/ZP;
			}
			realft(filteredprojection[l][i]-1,ZP,-1);    // IFFT
		}  // end of i-PROJECTION_Z
		
	 } // end of belta
    

	//step 3: back projection

	fp_h = (float*)malloc(X_SIZE * Y_SIZE * Z_SIZE * sizeof (float));
	for( i=0; i<X_SIZE; i++ )
		for( j=0; j<Y_SIZE; j++)
			for( k=0; k<Z_SIZE; k++)
			{
				fp_h[i*Y_SIZE*Z_SIZE+j*Z_SIZE+k] = filteredprojection[i][j][k];
			}
	 
	cudaMalloc((void**)&fp_d, X_SIZE * Y_SIZE * Z_SIZE * sizeof (float));
	cudaMemcpy(fp_d, fp_h, X_SIZE * Y_SIZE * Z_SIZE * sizeof (float), cudaMemcpyHostToDevice);
	
	cudaMalloc((void**)&CT_numbers_d, RECONSIZE_Z * RECONSIZE * RECONSIZE * sizeof (short int));
	cudaMemset((void*)CT_numbers_d, 0, RECONSIZE_Z * RECONSIZE * RECONSIZE * sizeof (short int));
	CT_numbers_h = (short int*)malloc(RECONSIZE_Z * RECONSIZE * RECONSIZE * sizeof (short int));
	
	dim3 dimBlock(THREADS, 1, 1);
	dim3 dimGrid((RECONSIZE_Z*RECONSIZE*RECONSIZE+THREADS-1)/THREADS, 1, 1);
    
	back_projection<<<dimGrid, dimBlock>>>(fp_d, CT_numbers_d, view_start, view_end, X_SIZE, Y_SIZE, Z_SIZE);

	cudaMemcpy(CT_numbers_h, CT_numbers_d, RECONSIZE_Z * RECONSIZE * RECONSIZE * sizeof (short int), cudaMemcpyDeviceToHost);
    
    gettimeofday(&tim, NULL);
    end = tim.tv_sec + (tim.tv_usec/1000000.0);
    
	fwrite(CT_numbers_h,sizeof(short int),RECONSIZE_Z * RECONSIZE * RECONSIZE,ptr_ct);
	
	fclose(ptr_proj);
    fclose(ptr_ct);
	free(fp_h);
	free(CT_numbers_h);
    

	printf("Execution time of FDK: %.2f seconds\n",end-begin);
} // end of main

/*********************************   FFT   ******************************************/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) 
	{
		if (j > i) 
		{
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) 
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) 
	{
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) 
		{
			for (i=m;i<=n;i+=istep) 
			{
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP

void realft(float data[], unsigned long n, int isign)
{
	void four1(float data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) 
	{
		c2 = -0.5;
		four1(data,n>>1,1);
	} 
	else 
	{
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) 
	{
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) 
	{
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} 
	else 
	{
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}

void cosft1(float y[], int n)
{
	void realft(float data[], unsigned long n, int isign);
	int j,n2;
	float sum,y1,y2;
	double theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;

	theta=PI/n;
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	sum=0.5*(y[1]-y[n+1]);
	y[1]=0.5*(y[1]+y[n+1]);
	n2=n+2;
	for (j=2;j<=(n>>1);j++) 
	{
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
		y1=0.5*(y[j]+y[n2-j]);
		y2=(y[j]-y[n2-j]);
		y[j]=y1-wi*y2;
		y[n2-j]=y1+wi*y2;
		sum += wr*y2;
	}
	realft(y,n,1);
	y[n+1]=y[2];
	y[2]=sum;
	for (j=4;j<=n;j+=2) 
	{
		sum += y[j];
		y[j]=sum;
	}
}

#undef PI
