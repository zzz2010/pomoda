import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import weka.core.Instances;


public class Domain2PWM {

	/**
	 * @param args
	 */
	
	static HashMap<Integer, Instances> loadTrainFile(String trainfile)
	{
		//load file for each AA, and do K-mean cluster for the discretize
		
		return null;	
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
			
			
		}
		catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "Domain2PWM", options );
			return;
		}
	}

}
