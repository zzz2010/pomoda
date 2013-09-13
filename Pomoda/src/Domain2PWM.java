import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.utils.ChangeVetoException;

import weka.clusterers.XMeans;
import weka.core.Attribute;
import weka.core.EuclideanDistance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.lazy.KStar;
import weka.classifiers.meta.AdaBoostM1;
import weka.classifiers.trees.J48graft;
import weka.classifiers.trees.RandomForest;

public class Domain2PWM {

	/**
	 * @param args
	 */
	static Instance Vec2Instance(String vec)
	{
		 double[] values = new double[4] ;
		 String[] comps=vec.split(",");
		 double sum=0;
		 for (int j = 0; j < 4; j++) {
			 values[j]=Double.parseDouble(comps[j]);
			 sum+=values[j];
		}	
		 //normalize to sum 1
		 for (int j = 0; j < 4; j++) {
			 values[j]=values[j]/sum;
		 }
		 Instance instance = new Instance(1, values);
		 
		 return instance;
		
	}
	static Instances Vecs2Instances(ArrayList<String> Vecs)
	{
		String ACGT="ACGT";
		FastVector attrList=new FastVector(4);
		for (int i = 0; i < 4; i++) {
			Attribute temp=new Attribute( ACGT.substring(i,i+1));
			attrList.addElement(temp);
		}
		Instances insts=new Instances("ACGT", attrList,Vecs.size());
		for (int i = 0; i < Vecs.size(); i++) {
			String vec = Vecs.get(i);
			Instance instance = Vec2Instance(vec);
			 insts.add(instance);
		}
		
		return insts;
	}
	
	static weka.clusterers.XMeans buildXMeanClusters(ArrayList<String> Vecs)
	{
		weka.clusterers.XMeans cluster=new XMeans();
		
		try {
			cluster.setMaxNumClusters(20);
			cluster.setMinNumClusters(4);
			Instances insts = Vecs2Instances(Vecs);
			cluster.buildClusterer(insts);
			
			Enumeration enumIt = insts.enumerateInstances();
			Instances centers = cluster.getClusterCenters();
			EuclideanDistance ED=new EuclideanDistance(insts);
			double sum=0;
			int count=0;
			while(enumIt.hasMoreElements())
			{
				Instance inst=(Instance) enumIt.nextElement();
				int classid=cluster.clusterInstance(inst);
				double dist = ED.distance(inst, centers.instance(classid))/(float)Math.sqrt(2);
				sum+=dist;
				count+=1;
			}
			System.out.println("Average Distance:"+sum/count);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return cluster;
	}
	
	static Attribute createWekaLabel_Vecs(weka.clusterers.XMeans cluster)
	{
		
		Instances insts = cluster.getClusterCenters();
		FastVector values = new FastVector();
		for (int i = 0; i < insts.numInstances(); i++) {
			values.addElement(insts.instance(i).toString());
		}
		return new Attribute("label",values);

	}
	
	static String AAstring="ABCDEFGHIKLMNPQRSTVWXYZ";
	
	static Attribute createWekaAttribute_AA(String name)
	{
		
		FastVector values = new FastVector();
		for (int i = 0; i < AAstring.length(); i++) {
			values.addElement(String.valueOf(AAstring.charAt(i)));
		}
		return new Attribute(name,values);
		
	}
	
	static Instances CreateInstances( ArrayList<String> AAs,  ArrayList<String> Vecs,int lastPWMColID,HashMap<String,Integer> AAmapping) throws Exception
	{
		  //do cluster for Vecs, parse instances object
		  FastVector attrList=new FastVector(AAs.get(0).length()+1);
		  for (int i = 0; i < AAs.get(0).length(); i++) {
			  attrList.addElement(createWekaAttribute_AA(String.valueOf(i+1)));
		}
		  //build cluster
		  weka.clusterers.XMeans cluster=buildXMeanClusters(Vecs);
		  Attribute label=createWekaLabel_Vecs(cluster);
		  System.out.println("PWM column "+lastPWMColID+"    cluster number:"+cluster.getClusterCenters().numInstances());
		  attrList.addElement(label);
		  Instances traindata=new Instances(String.valueOf(lastPWMColID), attrList,AAs.size());
		  //add instance to instances
		  for (int i = 0; i < AAs.size(); i++) {
			double[] values = new double[traindata.numAttributes()] ;
			String AAinst = AAs.get(i);
			for (int j = 0; j <AAinst.length(); j++) {
				values[j]=AAmapping.get(AAinst.substring(j, j+1));
			}
			Instance vecInst = Vec2Instance(Vecs.get(i));
			values[traindata.numAttributes()-1]=cluster.clusterInstance(vecInst);
			Instance instance = new Instance(1, values);
			traindata.add(instance);
		}
		  traindata.setClassIndex(traindata.numAttributes()-1);
		  return traindata;
	}
	
	static HashMap<String, HashMap<Integer,Instance>> loadTestFile(String testfile)
	{
		//prepare AA mapping table
		HashMap<String,Integer> AAmapping=new HashMap<String,Integer>();
		for (int i = 0; i < AAstring.length(); i++) {
			AAmapping.put(AAstring.substring(i,i+1), i);
		}
		//read the test file
		HashMap<String, HashMap<Integer,Instance>> testData=new HashMap<String, HashMap<Integer,Instance>>();
		   try {
	        	  BufferedReader sr = new BufferedReader(new FileReader(new File(testfile)));
	        	  String line = "";
	        	  while ((line = sr.readLine()) != null)
	    			{
	          		  String[] comps = line.trim().split("[ |\t]+");
	          		  if(comps.length<3)
	          			  continue;
	          		  String proteinName=comps[2];
	          		  int pwmColid=Integer.parseInt(comps[0]);
	          		  String AAinst=comps[1];
	          		double[] values = new double[AAinst.length()+1] ;
	    			for (int j = 0; j <AAinst.length(); j++) {
	    				values[j]=AAmapping.get(AAinst.substring(j, j+1));
	    			}
	    			values[AAinst.length()]=-1;
	    			Instance instance = new Instance(1, values);
	    			if(!testData.containsKey(proteinName))
	    				testData.put(proteinName, new HashMap<Integer,Instance>());
	    			testData.get(proteinName).put(pwmColid, instance);
	    			
	    			}
	        	  
		   }
	        catch(Exception e)
	        {
	        		e.printStackTrace();
	        }
		
		
		return testData;
	}
	
	static HashMap<Integer, Instances> loadTrainFile(String trainfile)
	{
		Random rrand=new Random();
		//prepare AA mapping table
		HashMap<String,Integer> AAmapping=new HashMap<String,Integer>();
		for (int i = 0; i < AAstring.length(); i++) {
			AAmapping.put(AAstring.substring(i,i+1), i);
		}
		HashMap<Integer, Instances> TrainDataCollection=new HashMap<Integer, Instances>();
		//load file for each AA, and do K-mean cluster for the discretize
        String line = "";
        try {
        	  BufferedReader sr = new BufferedReader(new FileReader(new File(trainfile)));
        	  ArrayList<String> AAs=new ArrayList<String>();
        	  ArrayList<String> Vecs=new ArrayList<String>();
        	  int  lastPWMColID=-1;
        	  int proteinId=0;
        	  while ((line = sr.readLine()) != null)
  			{
        		  String[] comps = line.trim().split("[ |\t]+");
        		  if(comps.length<3)
        			  continue;
        		  int PWMColID=Integer.parseInt(comps[0]);
        		  if(lastPWMColID!=PWMColID&&lastPWMColID!=-1)
        		  {
        			  Instances traindata = CreateInstances(AAs,Vecs,lastPWMColID,AAmapping);
        			  TrainDataCollection.put(lastPWMColID, traindata);
        			 
        			  //clear up AAs and Vecs
        			  proteinId=0;
        			  AAs.clear();
        			  Vecs.clear();
        		  }
        		  lastPWMColID=PWMColID;
        		  proteinId++;
        		  if(rrand.nextDouble()<0.1&&splitTrainTest) //10% test set
        		  {
        			  if(filterProtein.size()<proteinId*0.1)
        				  filterProtein.add(comps[3]);
        			  continue;
        		  }
        		  AAs.add(comps[1]);
        		  Vecs.add(comps[2]);
  			}
        	  //add last one 
        	  Instances traindata = CreateInstances(AAs,Vecs,lastPWMColID,AAmapping);
		  TrainDataCollection.put(lastPWMColID, traindata);
        	  
        }
        catch(Exception e)
        {
        		e.printStackTrace();
        }
		return TrainDataCollection;	
	}
	
	static HashSet<String> filterProtein=new HashSet<String>();
	static boolean splitTrainTest=true; //need the 4th column , can split the trainset and testset
	public static void main(String[] args) {
		
		// TODO Auto-generated method stub
		Options options = new Options();
		options.addOption("pwm", true, "PWM of combining all the training PWMs");
		options.addOption("train", true, "amino acid to PWM column file");
		options.addOption("test", true, "testing domain file");
		options.addOption("model", true, "training model object file");
		options.addOption("validate", true, "PWM of test protein sequences");
		
		CommandLineParser parser = new GnuParser();
		PWM FBP=null;
		List<PWM> validatePWM=null;
		String trainfile="";
		String testfile="";
		String modelfile="";
		try {
			CommandLine cmd = parser.parse( options, args);
			if(cmd.hasOption("pwm"))
			{
				FBP=common.LoadPWMFromFile(cmd.getOptionValue("pwm")).get(0);
			}
			if(cmd.hasOption("validate"))
			{
				validatePWM=common.LoadPWMFromFile(cmd.getOptionValue("validate"));
			}
			if(cmd.hasOption("train"))
			{
				trainfile=cmd.getOptionValue("train");
			}
			else
			{
				throw new ParseException("must provide the train file");
			}
			if(cmd.hasOption("model"))
			{
				modelfile=cmd.getOptionValue("model");
			}
			if(cmd.hasOption("test"))
			{
				testfile=cmd.getOptionValue("test");
			}
			

			
			HashMap<Integer, Instances> trainData=loadTrainFile(trainfile);
			HashMap<Integer, Classifier> trainModels=new HashMap<Integer, Classifier>();
			//build classifier for each PWM column
			System.out.println("===============Building Classification Models==================");
			for (Integer PWMcol : trainData.keySet()) {
				//Classifier Modeler=new weka.classifiers.lazy.KStar();
				Classifier Modeler=new RandomForest();
				//Classifier Modeler=new weka.classifiers.functions.MultilayerPerceptron();
				Instances data = trainData.get(PWMcol);
				Evaluation eval1 = new Evaluation(data);
				Modeler.buildClassifier(data);
				eval1.evaluateModel(Modeler, data);
				//eval1.crossValidateModel(Modeler, data, 20, new Random(1));//data.numInstances() leave one out
				System.out.println("PWM column "+PWMcol);
				System.out.println(eval1.toClassDetailsString());
				Modeler.buildClassifier(data);
				trainModels.put(PWMcol, Modeler);
			}
			HashMap<String, PWM> predictResults=new HashMap<String, PWM>();
			if(testfile!="")
			{
				HashMap<String, HashMap<Integer,Instance>> testData=loadTestFile(testfile);
				if(FBP!=null)
				{
					//output file
					BufferedWriter writer = new BufferedWriter(new FileWriter(testfile+".pwm"));
					//construct predicted PWM
					for (String proteinName : testData.keySet()) {
						PWM predictPWM=FBP.Clone();
						predictPWM.Name=proteinName;
						HashMap<Integer, Instance> PWMcolAAs = testData.get(proteinName);
						for (Integer PWMcol : PWMcolAAs.keySet()) {
							Classifier Modeler=trainModels.get(PWMcol);
							Instance testInst = PWMcolAAs.get(PWMcol);
							Instances dataSet = trainData.get(PWMcol);
							testInst.setDataset(dataSet);
							double classId=Modeler.classifyInstance(testInst);
							String predictVec=dataSet.classAttribute().value((int) classId);
							String[] toks=predictVec.split(",");
							for (int i = 0; i < 4; i++) {
								predictPWM.setWeight(PWMcol, i, Double.parseDouble(toks[i]));
							}
						}
						predictPWM=predictPWM.trim();
						predictPWM.Consensus(true);
						writer.write(predictPWM.toString()+"\n");
						predictResults.put(proteinName, predictPWM);
					}
					writer.close();
				}
				
			}
			
			if(validatePWM!=null)
			{
				double sumDiv=0;
				int count=0;
				int proteinId=0;
				System.out.println("===============Validation Result (PWM Divergence)==================");
				for (PWM pwm : validatePWM) {
					String name=pwm.Name.split("\\|")[0];
					proteinId++;
					if(!filterProtein.contains(name)&&splitTrainTest)
						continue;
					if(predictResults.containsKey(name))
					{
						PWM pwm2=predictResults.get(name);
						double pwmdiv=common.PWM_Divergence(pwm, pwm2);
						System.out.println(name+"\t"+pwmdiv);
						sumDiv+=pwmdiv;
						count+=1;
					}
				}
				if(count>0)
				System.out.println("Average PWM Divergence is "+sumDiv/count);
			}
			
			
		}
		catch (ParseException e) {
			// TODO Auto-generated catch block
			System.err.println(e.getMessage());
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "Domain2PWM", options );
			return;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
