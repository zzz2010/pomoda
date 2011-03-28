import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import org.apache.commons.lang.ArrayUtils;


import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;
import umontreal.iro.lecuyer.probdist.NormalDist;
import umontreal.iro.lecuyer.util.Num;

import edu.stanford.rsl.jpop.*;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;


public class NegBinFunction implements GradientOptimizableFunction {

	
	ArrayList<Integer> Rl;
	ArrayList<Double> Ez;
	double DampNegBin;
	int Num_Thread=1;
	
	public NegBinFunction(ArrayList<Integer> Rl_stat,ArrayList<Double> Ez_stat,double Damp)
	{
		Rl=Rl_stat;
		Ez=Ez_stat;
		DampNegBin=Damp;
	}

	public double[] run(double[] preNegBinpara)
	{
		
		double[] NegBinPar=preNegBinpara;
		if(NegBinPar==null)
			NegBinPar=initNegBinParamsMixture();
		int[] iflag=new int[] {0} ;
		
		//compute FG model first
		double[] NegBinPar1=new double[]{NegBinPar[2],NegBinPar[3]};
		FunctionOptimizer opt=new FunctionOptimizer(2);
		opt.setInitialX(NegBinPar1);
		opt.setOptimizationMode(OptimizationMode.Gradient);
		opt.setFunctionController(new SimpleFunctionController());
//		opt.setFunctionController(new ParallelFunctionController());
		opt.setFunctionAssembler(new AdditiveFunctionAssembler());
		opt.optimizeFunction(this);
		NegBinPar1=opt.getOptimum();
		
		
		//compute BG model
		for (int i = 0; i < Ez.size(); i++) {
			Ez.set(i, 1-Ez.get(i));
		}
		double[] NegBinPar0=new double[]{NegBinPar[0],NegBinPar[1]};
		opt=new FunctionOptimizer(2);
		opt.setInitialX(NegBinPar1);
		opt.setOptimizationMode(OptimizationMode.Gradient);
		opt.setFunctionController(new SimpleFunctionController());
//		opt.setFunctionController(new ParallelFunctionController());
		opt.setFunctionAssembler(new AdditiveFunctionAssembler());
		opt.optimizeFunction(this);
		NegBinPar0=opt.getOptimum();
		

		return new double[]{NegBinPar0[0],NegBinPar0[1],NegBinPar1[0],NegBinPar1[1]};
		
	}
	
	
	double [] initNegBinParamsMixture()
	{
		if(DampNegBin>0)
		{
			for (int i = 0; i < Ez.size(); i++) {
				Ez.set(i, Ez.get(i) * (1 - DampNegBin) + 0.5 * (DampNegBin));
			}
		}
		else
		{
			ArrayList<Integer> sortedRl=new ArrayList<Integer>(Rl);
			Collections.sort(sortedRl);
			for (int i = 0; i < Ez.size(); i++) {
				if(Rl.get(i)>sortedRl.get((int)(0.9*Rl.size())))
					Ez.set(i, 1.0);
				else
					Ez.set(i, 0.0);
				
			}
		}
		
		double[] NegBinPar=new double[4];
		double A0,A1,B0,B1;
		double myM=0;
		double myV=0;
		double sumREz=0;
		double sumEz=0;
		double sumR2Ez=0;
		
		double reduceFactor=1;
		ArrayList<Double> rl=new ArrayList<Double>();
		for (int i = 0; i < Rl.size(); i++) {
//			sumREz+=Rl.get(i)*Ez.get(i)/reduceFactor;
			sumEz+=Ez.get(i);
//			sumR2Ez+=Rl.get(i)*Rl.get(i)*Ez.get(i)/(reduceFactor*reduceFactor);
			if(Ez.get(i)>0.5)
				rl.add((double) Rl.get(i));
		}
		
		double[] ss=NormalDist.getMLE(ArrayUtils.toPrimitive(rl.toArray(new Double[1])),rl.size());
		myM=ss[0];//(sumREz/sumEz);
		myV=ss[1]*ss[1];//sumR2Ez/sumEz-myM*myM;
		B1=myM/myV;
		  if (B1 > 0.9) {
	            B1 =0.9;
	        }
		
		A1=myM*B1/(1-B1);
		
		
		sumREz=0;
		sumEz=0;
		sumR2Ez=0;
		rl.clear();
		for (int i = 0; i < Rl.size(); i++) {
//			sumREz+=Rl.get(i)*(1-Ez.get(i));
			sumEz+=(1-Ez.get(i));
//			sumR2Ez+=Rl.get(i)*Rl.get(i)*(1-Ez.get(i));
			if(Ez.get(i)<0.5)
				rl.add((double) Rl.get(i));
		}
		ss=NormalDist.getMLE(ArrayUtils.toPrimitive(rl.toArray(new Double[1])),rl.size());
		myM=ss[0];//(sumREz/sumEz);
		myV=ss[1]*ss[1];//sumR2Ez/sumEz-myM*myM;
		B0=myM/myV;
		  if (B0 > 0.9) {
	            B0 =0.9;
	        }
		A0=myM*B0/(1-B0);
		
		NegBinPar[1]=qlogis(B0);
		NegBinPar[0]=Math.log(A0);
		NegBinPar[3]=qlogis(B1);
		NegBinPar[2]=Math.log(A1);
		
		return NegBinPar;
	}
	
	public static double qlogis(double p)
	{
		return  Math.log(p/(1-p));
	}
	
	public static double plogis(double q)
	{
		double eq=Math.exp(q);
		return eq/(1+eq);
	}
	@Override
	public void setNumberOfProcessingBlocks(int number) {
		Num_Thread=number;
		
	}
	@Override
	public int getNumberOfProcessingBlocks() {
		// TODO Auto-generated method stub
		return Num_Thread;
	}
	@Override
	public double evaluate(double[] x, int block) {

		if(x[0]>200)
			x[0]=200;
		if(x[0]<-200)
			x[0]=-200;
//		if(x[1]>200)
//			x[1]=200;
//		if(x[1]<-200)
//			x[1]=-200;
		 double A = Math.exp(x[0]);
	      double  logitB = x[1];
	   // System.out.println("block:"+block);
	      double sum= 1e+300;
	        if ((A < 1e+200) && (A > common.DoubleMinNormal)) {
	        	 sum=0;
	   		  int start=(int) Math.ceil(block*Rl.size()/(double)getNumberOfProcessingBlocks());
			  int end=(int) Math.floor((block+1)*Rl.size()/(double)getNumberOfProcessingBlocks());
	    		for (int i = start; i < end; i++) {
	    			//System.out.println(Rl.get(i));
	    			sum+=-Ez.get(i)*(Num.lnGamma(A+Rl.get(i))-Num.lnGamma(A)+A*logitB-(A+Rl.get(i))*Math.log(1+Math.exp(logitB)));
	    		}
	           
	        }
	        //Num.lnGamma(A+Rl.get(i))-Num.lnGamma(A)
	        if(Double.isNaN(sum)||Double.isNaN(A)||Double.isInfinite(sum)||Double.isInfinite(A))
	        	sum= 1e+300;
	        	
	      return sum;
	}
	@Override
	public double[] gradient(double[] x, int block) {

		if(x[0]>200)
			x[0]=200;
		if(x[0]<-200)
			x[0]=-200;
//		if(x[1]>200)
//			x[1]=200;
//		if(x[1]<-200)
//			x[1]=-200;
		 double A = Math.exp(x[0]);
	      double  logitB = x[1];
		      double  B = plogis(logitB);
		  	double gLogA=0,gLogitB=0;
		  int start=(int) Math.ceil(block*Rl.size()/(double)getNumberOfProcessingBlocks());
		  int end=(int) Math.floor((block+1)*Rl.size()/(double)getNumberOfProcessingBlocks());
			for (int i = start; i < end; i++) {
				gLogA+=-A*(Ez.get(i)*(Num.digamma(Rl.get(i)+A)-Num.digamma(A)+logitB-Math.log(1 + Math.exp(logitB))));
				gLogitB+=-Ez.get(i)*(A-(A+Rl.get(i))*B);
			}
			
			//System.out.println(gLogA+"|||"+gLogitB);
		   return new double[]{gLogA,gLogitB};
	}
	
	
	public static void main(String[] args) {
		NegativeBinomialDist dist1=new NegativeBinomialDist(3000,0.9);	
		NegativeBinomialDist dist0=new NegativeBinomialDist(10,0.4);	
		Random r=new Random();
		double prior=0.7;
		ArrayList<Integer> Rl=new ArrayList<Integer>();
		for (int i = 0; i < 10000; i++) {
			if(r.nextDouble()<prior)
			{
				Rl.add(dist1.inverseFInt(r.nextDouble()));
			}
			else
			{
				Rl.add(dist0.inverseFInt(r.nextDouble()));
			}
		}
		ArrayList<Integer> sortedRl=new ArrayList<Integer>(Rl);
		Collections.sort(sortedRl);
		ArrayList<Double> Ez=new ArrayList<Double>();
		for (int i = 0; i < Rl.size(); i++) {
			if(Rl.get(i)>sortedRl.get((int)(Rl.size()*0.9)))
				Ez.add(1.0);
			else
				Ez.add(0.0);
				
		}
		double[] oldparas=null;
		while(true)
		{
			NegBinFunction solver=new NegBinFunction(Rl, Ez, 0);
			double Prior_Ez=0;
			for (int i = 0; i < Rl.size(); i++) {
				Prior_Ez+=Ez.get(i);
			}
			Prior_Ez/=Ez.size();
			double[] paras=solver.run(oldparas);
			NegativeBinomialDist DnaseBG=new NegativeBinomialDist(Math.exp( paras[0]),NegBinFunction.plogis( paras[1]));
			NegativeBinomialDist DnaseFG=new NegativeBinomialDist(Math.exp( paras[2]),NegBinFunction.plogis( paras[3]));
			for (int i = 0; i < Rl.size(); i++) {
				double r1=Prior_Ez*(DnaseFG.prob(Rl.get(i)))/((1-Prior_Ez)*(DnaseBG.prob(Rl.get(i))));
				double pl=r1/(1+r1);
				if(Double.isNaN(pl))
					pl=Prior_Ez;
				Ez.set(i, pl);
			}
			System.out.println(DnaseFG);
			System.out.println(DnaseBG);
			System.out.println(Prior_Ez);
			if(Arrays.equals(oldparas, paras))
				break;
			oldparas=paras;

		}
		
	}
	
	}


