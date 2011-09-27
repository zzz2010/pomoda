import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

import javax.imageio.ImageIO;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.bio.*;
import org.biojava.bio.dist.*;
import org.biojava.bio.dp.*;
import org.biojava.bio.gui.*;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.io.CrossProductTokenization;
import org.biojava.bio.seq.io.NameTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.bio.symbol.Symbol;


public class WMPanel extends JPanel {
	 private GapPWM wm;
     private DistributionLogo[] logos;
     HashMap<Integer,ArrayList<Integer>> groupList=new HashMap<Integer,ArrayList<Integer>>(10);
     HashMap<Integer,Double> groupBits=new HashMap<Integer, Double>();
     
     int width=40;
	 int height=120;
	 static String outputPrefix="";
	 
	 public static BufferedImage componentToBufferedImage(JComponent
			 component,
			 Dimension dim) throws IOException {
		 component.setDoubleBuffered(false);
		 component.setSize(dim);
		 component.addNotify();
			 BufferedImage backBuffer = new BufferedImage (dim.width, 
			 dim.height, BufferedImage.TYPE_INT_RGB);
			 Graphics g = backBuffer.createGraphics();
			 g.setColor(Color.white);
			 g.fillRect(0, 0, dim.width, dim.height);
			 g.setClip(0, 0, dim.width, dim.height);
			component.print(g);
			 
			g.dispose();
			 return backBuffer;

			 } // end componentToBufferedImage method
	 
	 
	public  void SaveLogoImage(String pngfile)
	{
		BufferedImage image =(BufferedImage) this.createImage(this.getPreferredSize().width, this.getPreferredSize().height);
//		BufferedImage image;
		try {
			image = componentToBufferedImage(this,this.getPreferredSize());
			
			ImageIO.write(image, "png", new File(pngfile ));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
	

	}
 	protected void paintComponent(Graphics g) {
 		 
 		Graphics2D g2 = (Graphics2D) g;
			
 		int groupsize=Math.max(1,groupList.size());//when there is no dependgroup set to 1
//      sorted by information bit
 		ArrayList<Double> sortedBit=null;
 		if(groupsize>1)
 		{
 			sortedBit=new ArrayList<Double>(groupBits.values());
 			Collections.sort(sortedBit);
 			
 		}
 		
 		
 		int gheight=(int) ((double)height/groupsize*0.8);//make some room for label
		int gid=groupList.size();
		for(Entry<Integer, ArrayList<Integer>> pair:groupList.entrySet())
		{
			ArrayList<Integer> poslist=pair.getValue();
			int lastx=-1;
			if(sortedBit!=null)
			gid=Collections.binarySearch(sortedBit, groupBits.get(pair.getKey()))+1;
			for (int i = 0; i < poslist.size(); i++) {
				//draw vertical line
				int x=(int) ((poslist.get(i)+0.5)*width);
				 Stroke drawingStroke = new BasicStroke(3, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{9}, 0);
				 g2.setStroke(drawingStroke);
				g2.drawLine ( x, height, x, height+gid*gheight );
				if(i==poslist.size()/3)
				{
					JLabel label=new JLabel(groupBits.get(pair.getKey()).toString().substring(0, Math.min(5,groupBits.get(pair.getKey()).toString().length())));
					label.setBounds(x, height+gid*gheight, width, height/10);
					add(label);
				}
				if(lastx>0)
				{
					g2.setStroke(new BasicStroke(2.0f));
					g2.drawLine ( lastx, height+gid*gheight, x, height+gid*gheight );
				}
				lastx=x;
			}
		
		}
		
 
	}
     public WMPanel(GapPWM wm) {
         super();
         
         this.wm = wm;
         setBackground(Color.white);
         groupList.clear();
         groupBits.clear();
         RenderingHints hints = new 
RenderingHints(RenderingHints.KEY_ANTIALIASING, 
RenderingHints.VALUE_ANTIALIAS_ON);

    	 char[] ACGT=new char[]{'A','C','G','T'};
         try {
             //setLayout(new GridLayout(1, wm.columns()));
        	 setLayout(null);
        	 this.addNotify();
        	 this.setDoubleBuffered(false);
        	 HashSet<Integer> visited=new HashSet<Integer>();
             logos = new DistributionLogo[wm.columns()];
             for (int pos = 0; pos < wm.columns(); ++pos) {
                 Distribution dist = null;
                 int nwidth=1;
                 if(wm.GroupId[pos]==0)
                 {
                	 dist=wm.getColumn(pos);
                 }
                 else
                 {	
                	 if(groupList.containsKey(wm.GroupId[pos]))
                		 groupList.get(wm.GroupId[pos]).add(pos);
                	 else
                	 {
                		 groupList.put(wm.GroupId[pos], new ArrayList<Integer>());
                		 groupList.get(wm.GroupId[pos]).add(pos);
                	 }
                	 if(visited.contains(wm.GroupId[pos]))
                		 continue;
                	 
                	 HashMap<String,Double> dprobs=wm.Dgroup_DmerProb.get(wm.GroupId[pos]);
                	 int groupsize=dprobs.keySet().iterator().next().length();
                	 String[] Nstrs=new  String[groupsize-1];
                	 String temp="";
                	 int gid=0;
                	 for (int i = pos+1; i < wm.columns(); i++) {
                		 if(wm.GroupId[pos]!=wm.GroupId[i])
                			 temp+=" ";
                		 else
                		 {
                			 Nstrs[gid]=temp;
                			 temp="";
                			 gid++;
                		 }
					}
                	  
                	  SimpleAlphabet  groupAlphabet=new SimpleAlphabet();
                	  groupAlphabet.setName("joinDNA");
                	  
                	  for (int i = 0; i < Math.pow(4, groupsize); i++) {	
                		   temp="";
                		   int aid=i;
                		  for (int j = 0; j < groupsize; j++) {
							temp+=ACGT[aid%4];
							aid/=4;
							if(j<groupsize-1)
							temp+=Nstrs[j];
						}
                		  Symbol sym = AlphabetManager.createSymbol(temp,Annotation.EMPTY_ANNOTATION);
                		  groupAlphabet.addSymbol(sym);
                	  }
                	   groupAlphabet.putTokenization("token", new NameTokenization(groupAlphabet));
                	  dist = DistributionFactory.DEFAULT.createDistribution(groupAlphabet);

                	  Iterator<Symbol> iter=groupAlphabet.iterator();
                	  double sum=0;
                	  while(iter.hasNext())
                	  {
                		  Symbol sym=iter.next();
                		  double weight=dprobs.get("N");
                		  nwidth=sym.getName().length();
                		String key=sym.getName().replace(" ", "");  
                		if(dprobs.containsKey(key))
                			weight=dprobs.get(key);
                		sum+=weight;
                		dist.setWeight(sym, weight);
					}
                	  
                	  groupBits.put(wm.GroupId[pos], DistributionTools.bitsOfInformation(dist)) ;
                	  visited.add(wm.GroupId[pos]);
                 }
                
                 DistributionLogo dl = new DistributionLogo();
                 if(nwidth>1)
                	 dl.setScaleByInformation(false);
                 dl.setRenderingHints(hints);
                 dl.setBackground(Color.white);
                 
                 dl.setOpaque(false);
                 dl.setDistribution(dist);
                 dl.setPreferredSize(new Dimension(width*nwidth, height));
                 dl.setLogoPainter(new gapTextLogoPainter());
                 dl.setStyle(new DNAStyle());
                 dl.setBounds(pos*width, 0, width*nwidth, height);
                 dl.addNotify();
                 add(dl);
                 logos[pos] = dl;
             }
             this.setPreferredSize(new Dimension( width*wm.columns(), height*2));
            
         } catch (BioException ex) {
             throw new BioError(ex);
         }
     }

     public static void wmViewer(GapPWM wm, String message)
     {
    	 
         WMPanel wmv = new WMPanel(wm);
         if(java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadlessInstance())
         {
        	 wmv.SaveLogoImage(outputPrefix+wm.Name+ ".png");
        	 
         }
         else
         {
         JFrame frame = new JFrame("Weight matrix viewer" + ((message == 
null) ? "" : (" (" + message + ")")));
         frame.getContentPane().add(wmv);
        
         frame.pack();
         
         frame.setVisible(true);
         }
        
         
     }
     public static void main(String[] args) {
    	 int topN=100;
    	 String inputPWMfile="";
    	
    		Options options = new Options();
    		options.addOption("i", true, "input dependency PWM file");
    		options.addOption("N", true, "the number of PWMs used to generate png files");
    		options.addOption("prefix", true, "output directory");
    		String inputPWM;
    		int simlen=-1;
    		CommandLineParser parser = new GnuParser();
    		
    		try {
    			CommandLine cmd = parser.parse( options, args);
    			if(cmd.hasOption("i"))
    			{
    				inputPWMfile=cmd.getOptionValue("i");
    			}
    			else
    			{
    				throw new ParseException("no input fasta file");
    			}


    			if(cmd.hasOption("N"))
    			{
    				topN=Integer.parseInt(cmd.getOptionValue("N"));
    			}
 
    			if(cmd.hasOption("prefix"))
    			{
    				outputPrefix=cmd.getOptionValue("prefix");
    			}
    
    		} catch (ParseException e) {
    			// TODO Auto-generated catch block
    			HelpFormatter formatter = new HelpFormatter();
    			formatter.printHelp( "FastaMask", options );
    			return;
    		}
    	 LinkedList<PWM> pwms=common.LoadPWMFromFile(inputPWMfile);
    	 for (int i = 0; i < topN&&i<pwms.size(); i++) {
    		 if(pwms.get(i) instanceof GapPWM)
    		 wmViewer((GapPWM)pwms.get(i), "Dependency Motif Logo");
		}
    	
    	 
    	 
     }
     
}
