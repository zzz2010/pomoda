import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.utils.ChangeVetoException;

import weka.clusterers.XMeans;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;


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
			cluster.setMaxKMeans(20);
			cluster.setMinNumClusters(2);
			cluster.buildClusterer(Vecs2Instances(Vecs));
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
			values.addElement(AAstring.charAt(i));
		}
		return new Attribute(name,values);
		
	}
	
	static HashMap<Integer, Instances> loadTrainFile(String trainfile)
	{
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
        	  while ((line = sr.readLine()) != null)
  			{
        		  String[] comps = line.trim().split("[ |\t]+");
        		  if(comps.length<3)
        			  continue;
        		  int PWMColID=Integer.parseInt(comps[0]);
        		  if(lastPWMColID!=PWMColID&&lastPWMColID!=-1)
        		  {
        			  //do cluster for Vecs, parse instances object
        			  FastVector attrList=new FastVector(AAs.get(0).length()+1);
        			  for (int i = 0; i < AAs.get(0).length(); i++) {
        				  attrList.addElement(createWekaAttribute_AA(String.valueOf(i+1)));
					}
        			  //build cluster
        			  weka.clusterers.XMeans cluster=buildXMeanClusters(Vecs);
        			  Attribute label=createWekaLabel_Vecs(cluster);
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
        			  TrainDataCollection.put(lastPWMColID, traindata);
        			 
        			  //clear up AAs and Vecs
        			  AAs.clear();
        			  Vecs.clear();
        		  }
        		  lastPWMColID=PWMColID;
        		  AAs.add(comps[1]);
        		  Vecs.add(comps[2]);
  			}
        }
        catch(Exception e)
        {
        		e.printStackTrace();
        }
		return TrainDataCollection;	
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Options options = new Options();
		options.addOption("FBP", true, "PWM of combining all the training PWMs");
		options.addOption("train", true, "amino acid to PWM column file");
		options.addOption("test", true, "testing domain file");
		options.addOption("model", true, "training model object file");
		
		CommandLineParser parser = new GnuParser();
		PWM FBP=null;
		String trainfile="";
		String testfile="";
		String modelfile="";
		try {
			CommandLine cmd = parser.parse( options, args);
			if(cmd.hasOption("pwm"))
			{
				FBP=common.LoadPWMFromFile(cmd.getOptionValue("pwm")).get(0);
			}
			if(cmd.hasOption("train"))
			{
				trainfile=cmd.getOptionValue("train");
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
			
			
		}
		catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "Domain2PWM", options );
			return;
		}
	}

}
